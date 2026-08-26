#!/bin/bash
set -eo pipefail

# ==================== 配置区（根据你的实际路径修改） ====================
workdir=
INPUT_VCF=${workdir}
INCLUDE_BED=${workdir}/out_genmap/high_mappability.merge.bed
OUTDIR=${workdir}/
OUT_PREFIX=SNP.filter

# 工具路径
VCFTOOLS=
BGZIP=
TABIX=

mkdir -p ${OUTDIR}

# 公共过滤参数
COMMON="--minQ 30 --min-meanDP 5 --max-meanDP 150 --bed ${INCLUDE_BED} --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all"
# ==================== 梯度1 ====================
${VCFTOOLS} --gzvcf ${INPUT_VCF} ${COMMON} \
    --max-missing 0.9 --mac 1 \
    --stdout | ${BGZIP} -c > ${OUTDIR}/${OUT_PREFIX}.first.vcf.gz
${TABIX} -p vcf ${OUTDIR}/${OUT_PREFIX}.first.vcf.gz
# ==================== 梯度2 ====================
${VCFTOOLS} --gzvcf ${INPUT_VCF} ${COMMON} \
    --max-missing 0.5 --mac 1 \
    --stdout | ${BGZIP} -c > ${OUTDIR}/${OUT_PREFIX}.second.vcf.gz
${TABIX} -p vcf ${OUTDIR}/${OUT_PREFIX}.secod.vcf.gz

# ==================== 梯度3 ====================
${VCFTOOLS} --gzvcf ${INPUT_VCF} ${COMMON} \
    --max-missing 0.2 --mac 1 \
    --stdout | ${BGZIP} -c > ${OUTDIR}/${OUT_PREFIX}.third.vcf.gz
${TABIX} -p vcf ${OUTDIR}/${OUT_PREFIX}.third.vcf.gz
# ==================== 梯度4 ====================
${VCFTOOLS} --gzvcf ${INPUT_VCF} ${COMMON} \
    --max-missing 0.9 --mac 0.5 \
    --stdout | ${BGZIP} -c > ${OUTDIR}/${OUT_PREFIX}.fourth.vcf.gz
${TABIX} -p vcf ${OUTDIR}/${OUT_PREFIX}.fourth.vcf.gz

echo "All done. Output in ${OUTDIR}"
