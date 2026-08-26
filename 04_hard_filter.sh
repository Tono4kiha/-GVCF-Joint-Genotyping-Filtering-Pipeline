# ==================== 软件路径 ====================
gatk=
samtools=
bgzip=
tabix=

# ==================== 文件路径 ====================
refgenome=
vcf=
out=
tmp=
species=

# Java 内存选项
java_opts="-Xmx20g -Djava.io.tmpdir=${tmp}"

# ==================== 创建必要目录 ====================
mkdir -p ${out} ${tmp}

# ==================== 准备参考基因组索引 ====================
if [ ! -f ${refgenome}.fai ]; then
    echo "Creating .fai index for reference..."
    ${samtools} faidx ${refgenome}
fi
if [ ! -f ${refgenome%.fa}.dict ]; then
    echo "Creating .dict sequence dictionary for reference..."
    ${gatk} CreateSequenceDictionary -R ${refgenome} -O ${refgenome%.fa}.dict
fi

if [ ! -f ${vcf} ]; then
    echo "ERROR: Input VCF missing: ${vcf}"
    exit 1
fi

# ==================== 1. 提取 SNP 和 INDEL ====================
echo "===== Step 1: Extract SNPs and INDELs ====="
${gatk} --java-options "${java_opts}" SelectVariants \
    -R ${refgenome} -V ${vcf} -select-type SNP -O ${out}/${species}.snp.raw.vcf.gz

${gatk} --java-options "${java_opts}" SelectVariants \
    -R ${refgenome} -V ${vcf} -select-type INDEL -O ${out}/${species}.indel.raw.vcf.gz

wait

# ==================== 2. 硬过滤（打标签，不删除位点） ====================
echo "===== Step 2: Apply hard filters to SNPs ====="
${gatk} --java-options "${java_opts}" VariantFiltration \
    -R ${refgenome} \
    --variant ${out}/${species}.snp.raw.vcf.gz \
    -O ${out}/${species}.snp.filtered.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name QD2 \
    --filter-expression "QUAL < 30.0" --filter-name QUAL30 \
    --filter-expression "FS > 20.0" --filter-name FS20 \
    --filter-expression "MQ < 20.0" --filter-name MQ20 \
    --filter-expression "MQRankSum < -0.5" --filter-name MQRankSum-0.5 \
    --filter-expression "ReadPosRankSum < -0.5" --filter-name ReadPosRankSum-0.5

echo "===== Step 3: Apply hard filters to INDELs ====="
${gatk} --java-options "${java_opts}" VariantFiltration \
    -R ${refgenome} \
    --variant ${out}/${species}.indel.raw.vcf.gz \
    -O ${out}/${species}.indel.filtered.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name QD2 \
    --filter-expression "QUAL < 30.0" --filter-name QUAL30 \
    --filter-expression "FS > 20.0" --filter-name FS20 \
    --filter-expression "MQ < 20.0" --filter-name MQ20

wait

# ==================== 3. 提取 PASS 位点（zcat + awk + bgzip） ====================
echo "===== Step 4: Extract PASS variants via zcat + awk + bgzip ====="
zcat ${out}/${species}.snp.filtered.vcf.gz | awk '$1~/^#/ || $7=="PASS"' | ${bgzip} > ${out}/${species}.snp.final.hardfilter.vcf.gz
zcat ${out}/${species}.indel.filtered.vcf.gz | awk '$1~/^#/ || $7=="PASS"' | ${bgzip} > ${out}/${species}.indel.final.hardfilter.vcf.gz

wait

# ==================== 4. 用 tabix 建立索引（假设文件已排序） ====================
echo "===== Step 5: Index final VCFs with tabix ====="
${tabix} -p vcf ${out}/${species}.snp.final.hardfilter.vcf.gz
${tabix} -p vcf ${out}/${species}.indel.final.hardfilter.vcf.gz

echo "===== All done! ====="
ls -lh ${out}/${species}.*.final.hardfilter.vcf.gz*
