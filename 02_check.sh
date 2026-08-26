#!/bin/bash
# ============================================================
# 脚本名称: extract_qc_metrics.sh
# 功能: 从 VCF 提取硬过滤六指标 + 统计信息 + 处理历史
# 所有输出均写入文件，不截断，不预览
# ============================================================

# ==================== 配置区（请按需修改） ====================
VCF=""
OUTDIR="./check_results"   # 结果输出目录
zgrep=
bcftools=
# 创建输出目录
mkdir -p "${OUTDIR}"

# 定义各输出文件路径
GATK_HIST="${OUTDIR}/gatk_history.txt"
STATS_OUT="${OUTDIR}/bcftools_stats.txt"
SAMPLE_LIST="${OUTDIR}/sample_list.txt"
METRICS_FILE="${OUTDIR}/snp_QC_metrics.txt"
SUMMARY_FILE="${OUTDIR}/metrics_summary.txt"
LOG_FILE="${OUTDIR}/run.log"

# 所有屏幕输出同时写入总日志
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "===== QC Metrics Extraction Started at $(date) ====="
echo "Input VCF: $VCF"
echo "Output folder: $OUTDIR"
echo ""

# 检查 VCF 是否存在
if [ ! -f "$VCF" ]; then
    echo "ERROR: VCF file not found: $VCF"
    exit 1
fi

# --------------------------- 1. GATK 处理历史 ---------------------------
#echo ">>> Extracting GATK processing history..."
${zgrep} "^##GATKCommandLine" "$VCF" > "${GATK_HIST}" 2>/dev/null
#if [ -s "${GATK_HIST}" ]; then
    echo "    Saved to: ${GATK_HIST}"
#else
    echo "    No GATKCommandLine lines found; file will be empty."
#fi

# --------------------------- 2. bcftools stats 完整统计 ---------------------------
echo ">>> Running bcftools stats (full output)..."
${bcftools} stats "$VCF" > "${STATS_OUT}" 2>&1
echo "    Saved to: ${STATS_OUT}"

# --------------------------- 3. 样本列表 ---------------------------
echo ">>> Extracting sample list..."
${bcftools} query -l "$VCF" > "${SAMPLE_LIST}"
TOTAL_SAMPLES=$(wc -l < "${SAMPLE_LIST}")
echo "    Total samples: ${TOTAL_SAMPLES}"
echo "    Saved to: ${SAMPLE_LIST}"

# --------------------------- 4. 提取六个硬过滤指标 ---------------------------
echo ">>> Extracting QC metrics: QD, FS, SOR, MQ, MQRankSum, ReadPosRankSum..."
${bcftools} query -f '%QD\t%FS\t%SOR\t%MQ\t%MQRankSum\t%ReadPosRankSum\n' "$VCF" > "${METRICS_FILE}"
TOTAL_RECORDS=$(wc -l < "${METRICS_FILE}")
echo "    Total records: ${TOTAL_RECORDS}"
echo "    Saved to: ${METRICS_FILE}"

# --------------------------- 5. 各指标汇总统计 ---------------------------
echo ">>> Calculating summary statistics for each metric..."
awk -F'\t' '
BEGIN {
    fields[1]="QD"; fields[2]="FS"; fields[3]="SOR";
    fields[4]="MQ"; fields[5]="MQRankSum"; fields[6]="ReadPosRankSum";
}
{
    for(i=1;i<=6;i++) {
        if($i != "." && $i != "NaN") {
            count[i]++; sum[i]+=$i; sumsq[i]+=$i*$i;
            if($i < min[i] || count[i]==1) min[i]=$i;
            if($i > max[i] || count[i]==1) max[i]=$i;
        }
    }
}
END {
    for(i=1;i<=6;i++) {
        if(count[i] > 0) {
            mean = sum[i]/count[i];
            std = sqrt(sumsq[i]/count[i] - mean*mean);
            printf "%s:\tNon-missing=%d\tMean=%.3f\tStd=%.3f\tMin=%.3f\tMax=%.3f\n",
                   fields[i], count[i], mean, std, min[i], max[i];
        } else {
            printf "%s:\tNon-missing=0\n", fields[i];
        }
    }
}' "${METRICS_FILE}" > "${SUMMARY_FILE}"
echo "    Saved to: ${SUMMARY_FILE}"

# --------------------------- 完成 ---------------------------
echo ""
echo "===== QC Metrics Extraction Finished at $(date) ====="
echo "All outputs saved in: ${OUTDIR}"
echo "Log file: ${LOG_FILE}"
