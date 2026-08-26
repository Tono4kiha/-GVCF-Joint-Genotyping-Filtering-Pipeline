 #!/usr/bin/env Rscript
# ============================================================
# 功能: 读取 VCF QC 指标文件，绘制密度分布图并输出 PNG
# 用法:
#   Rscript plot_QC_metrics.R [input.txt] [output.png]
# 默认:
#   input.txt = ./snp_QC_metrics.txt
#   output.png = ./QC_distributions.png
# ============================================================

# 解析命令行参数
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {
  input_file <- args[1]
} else {
  input_file <- "snp_QC_metrics.txt"
}

if (length(args) >= 2) {
  output_png <- args[2]
} else {
  output_png <- "QC_distributions.png"
}

# 检查输入文件是否存在
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# 读取数据（制表符分隔，缺失值用 '.' 表示）
df <- read.table(input_file, header = FALSE, sep = "\t", na.strings = ".",
                 col.names = c("QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum"))

# 设置绘图设备
png(output_png, width = 1200, height = 800, res = 120)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

# ---- 辅助绘图函数 ----
plot_density <- function(x, main, xlab, log10_transform = FALSE, vline = NULL, vline_label = NULL) {
  # 去除 NA 和 Inf
  x_clean <- x[is.finite(x)]
  if (length(x_clean) == 0) {
    plot(1, type = "n", main = main, xlab = xlab, ylab = "Density")
    text(1, 1, "No valid data", col = "red")
    return()
  }
  
  if (log10_transform) {
    # 对正值取 log10，避免 log(0)
    x_clean <- x_clean[x_clean > 0]
    if (length(x_clean) == 0) {
      plot(1, type = "n", main = main, xlab = xlab, ylab = "Density")
      text(1, 1, "No positive data for log10", col = "red")
      return()
    }
    x_plot <- log10(x_clean)
    xlab <- paste("log10(", xlab, ")", sep = "")
  } else {
    x_plot <- x_clean
  }
  
  dens <- density(x_plot)
  plot(dens, main = main, xlab = xlab, ylab = "Density", col = "steelblue", lwd = 2)
  
  # 如果指定了参考线，绘制虚线
  if (!is.null(vline)) {
    if (log10_transform) {
      vline_val <- log10(vline)
    } else {
      vline_val <- vline
    }
    abline(v = vline_val, col = "darkred", lty = 2, lwd = 1.5)
    if (!is.null(vline_label)) {
      text(vline_val, max(dens$y) * 0.9, labels = vline_label, col = "darkred", pos = 4, cex = 0.8)
    }
  }
}

# ---- 绘制六个指标 ----
# 1. QD (低值过滤，人类阈值 2.0)
plot_density(df$QD, "QD (Quality by Depth)", "QD", vline = 2.0, vline_label = "Human: 2.0")

# 2. FS (高值过滤，通常用 log10 变换，人类阈值 60.0)
plot_density(df$FS, "FS (Fisher Strand)", "FS", log10_transform = TRUE,
             vline = 60.0, vline_label = "Human: 60.0")

# 3. SOR (高值过滤，log10 变换，人类阈值 3.0)
plot_density(df$SOR, "SOR (Strand Odds Ratio)", "SOR", log10_transform = TRUE,
             vline = 3.0, vline_label = "Human: 3.0")

# 4. MQ (低值过滤，人类阈值 40.0)
plot_density(df$MQ, "MQ (RMS Mapping Quality)", "RMS Mapping Quality",
             vline = 40.0, vline_label = "Human: 40.0")

# 5. MQRankSum (低值过滤，人类阈值 -12.5)
plot_density(df$MQRankSum, "MQRankSum", "MQRankSum",
             vline = -12.5, vline_label = "Human: -12.5")

# 6. ReadPosRankSum (低值过滤，人类阈值 -8.0)
plot_density(df$ReadPosRankSum, "ReadPosRankSum", "ReadPosRankSum",
             vline = -8.0, vline_label = "Human: -8.0")

dev.off()

cat("Distribution plot saved to: ", output_png, "\n")
