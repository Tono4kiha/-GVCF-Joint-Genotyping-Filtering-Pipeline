#!/bin/bash
# ============================================================
# 按染色体拆分"合并 GVCF + 联合分型"脚本生成器
# 对应 combaine.sh 的两步流程（CombineGVCFs + GenotypeGVCFs），
# 本脚本按染色体逐一拆分：
#   一条染色体 -> 一个 interval 文件 -> 一个 qsub 子脚本
# 每个子脚本对指定染色体执行 CombineGVCFs + GenotypeGVCFs，
# 最后用 GatherVcfs 合并各染色体的 VCF 并建索引。
# 用法：bash generate_combine_per_chr.sh
# ============================================================

set -euo pipefail

# -------------------- 路径配置（沿用 combaine.sh）--------------------
workdir=
ref=
args_file=${workdir}/gvcf_args.txt
TMP=${workdir}/01_work/tmp
GATK=

outdir=${workdir}/01_work
per_chr_dir=${outdir}/per_chr
scriptdir=${outdir}/scripts
logdir=${outdir}/logs
intervaldir=${outdir}/intervals
prefix=

# -------------------- 准备工作 --------------------
mkdir -p ${TMP} ${per_chr_dir} ${scriptdir} ${logdir} ${intervaldir}

# 确保参考基因组及索引存在
if [ ! -f "${ref}" ]; then
    echo "错误: 未找到参考基因组 ${ref}"
    exit 1
fi
if [ ! -f "${ref}.fai" ]; then
    echo "错误: 未找到 ${ref}.fai，请先运行 samtools faidx ${ref}"
    exit 1
fi
if [ ! -f "${args_file}" ]; then
    echo "错误: 未找到样本列表 ${args_file}"
    exit 1
fi

# 获取染色体列表及对应的安全文件名（替换非字母数字为下划线，避免特殊字符）
chr_list=()
safe_chr_list=()
while IFS= read -r chr; do
    if [ -z "${chr}" ]; then
        continue
    fi
    chr_list+=("${chr}")
    safe_chr_list+=("${chr//[^A-Za-z0-9]/_}")
done < <(cut -f1 "${ref}.fai")
total_chr=${#chr_list[@]}

if [ "${total_chr}" -lt 1 ]; then
    echo "错误: 未从 ${ref}.fai 读取到任何染色体"
    exit 1
fi

echo "========================================"
echo "  参考基因组: ${ref}"
echo "  样本列表:   ${args_file}"
echo "  染色体总数: ${total_chr}"
echo "  将生成 ${total_chr} 个独立 qsub 脚本（每条染色体一个）"
echo "========================================"

# -------------------- 逐条染色体生成子脚本 --------------------
for ((i=0; i<total_chr; i++)); do
    chr=${chr_list[$i]}
    safe_chr=${safe_chr_list[$i]}

    # 生成该染色体的 interval 文件
    interval_file=${intervaldir}/${safe_chr}.list
    echo "${chr}" > ${interval_file}

    # 生成该染色体的 qsub 子脚本
    script_file=${scriptdir}/combine_chr_${safe_chr}.sh
    cat > ${script_file} <<EOF
#!/bin/bash

echo "===== 染色体: ${chr} ====="
echo "开始时间: \$(date)"

# 1) 合并该染色体的所有样本 GVCF
${GATK} --java-options "-Xmx20g -Djava.io.tmpdir=${TMP}" \\
    CombineGVCFs \\
    -R ${ref} \\
    -L ${interval_file} \\
    --arguments_file ${args_file} \\
    -O ${per_chr_dir}/${prefix}_merge.${safe_chr}.g.vcf.gz

if [ \$? -ne 0 ]; then
    echo "染色体 ${chr} CombineGVCFs 失败！" >&2
    exit 1
fi

# 2) 联合基因分型
${GATK} --java-options "-Xmx20g -Djava.io.tmpdir=${TMP}" \\
    GenotypeGVCFs \\
    -R ${ref} \\
    -V ${per_chr_dir}/${prefix}_merge.${safe_chr}.g.vcf.gz \\
    -O ${per_chr_dir}/${prefix}.${safe_chr}.vcf.gz

if [ \$? -ne 0 ]; then
    echo "染色体 ${chr} GenotypeGVCFs 失败！" >&2
    exit 1
fi

echo "染色体 ${chr} 完成"
echo "结束时间: \$(date)"
EOF

    chmod +x ${script_file}
    echo "  [染色体 ${chr}] -> ${script_file}"
done

# -------------------- 生成合并脚本 --------------------
merge_script=${scriptdir}/merge_chr.sh
{
    echo "#!/bin/bash"
    echo ""
    echo "echo \"开始合并 ${total_chr} 条染色体的 VCF\""
    echo "echo \"开始时间: \$(date)\""
    echo ""
    echo "${GATK} --java-options \"-Xmx20g -Djava.io.tmpdir=${TMP}\" \\"
    echo "    GatherVcfs \\"
    for ((i=0; i<total_chr; i++)); do
        echo "    -I ${per_chr_dir}/${prefix}.${safe_chr_list[$i]}.vcf.gz \\"
    done
    echo "    -O ${outdir}/${prefix}.vcf.gz"
    echo ""
    echo "# 建索引"
    echo "${GATK} IndexFeatureFile -I ${outdir}/${prefix}.vcf.gz"
    echo ""
    echo "echo \"合并完成: \$(date)\""
} > ${merge_script}

chmod +x ${merge_script}

# -------------------- 生成一键提交脚本 --------------------
submit_all=${scriptdir}/submit_all.sh
{
    echo "#!/bin/bash"
    echo "# 一键提交所有染色体 CombineGVCFs+GenotypeGVCFs 任务"
    echo "echo \"提交 ${total_chr} 个染色体任务...\""
    echo ""
    for ((i=0; i<total_chr; i++)); do
        echo "qsub ${scriptdir}/combine_chr_${safe_chr_list[$i]}.sh"
        echo "echo \"  已提交: ${chr_list[$i]}\""
    done
    echo ""
    echo "echo \"全部提交完毕！\""
    echo "echo \"所有任务完成后，请手动提交合并任务:\""
    echo "echo \"  qsub ${merge_script}\""
} > ${submit_all}

chmod +x ${submit_all}

# -------------------- 输出汇总 --------------------
echo ""
echo "========================================"
echo "  Dry Run 完成！生成文件汇总："
echo "========================================"
echo ""
echo "  interval 文件目录: ${intervaldir}/"
echo "  任务脚本目录:      ${scriptdir}/"
echo "  日志输出目录:      ${logdir}/"
echo ""
echo "  染色体任务脚本 (${total_chr} 个):"
for ((i=0; i<total_chr; i++)); do
    echo "    ${scriptdir}/combine_chr_${safe_chr_list[$i]}.sh"
done
echo ""
echo "  合并脚本:"
echo "    ${merge_script}"
echo ""
echo "  一键提交脚本:"
echo "    ${submit_all}"
echo ""
echo "========================================"
echo "  使用方式："
echo "========================================"
echo ""
echo "  方式一：一键提交所有染色体任务"
echo "    bash ${submit_all}"
echo ""
echo "  方式二：手动逐个提交"
for ((i=0; i<total_chr; i++)); do
    echo "    qsub ${scriptdir}/combine_chr_${safe_chr_list[$i]}.sh"
done
echo ""
echo "  全部完成后提交合并任务："
echo "    qsub ${merge_script}"
echo "========================================"
