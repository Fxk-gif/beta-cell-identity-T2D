library(ggplot2)

base <- "C:/Users/45902/WorkBuddy/2026-09-20-15-12-44/r3_inputs"
outdir <- "C:/Users/45902/WorkBuddy/2026-09-20-15-12-44/figures"

# 读全基因 logFC
fc <- read.csv(file.path(base, "volcano_input.csv"))
names(fc) <- c("gene", "logFC", "padj")
fc <- fc[order(-fc$logFC), ]  # 降序

# 读靶基因集
sets <- read.csv(file.path(base, "tf_target_sets.csv"))
set_names <- unique(sets$set)

# GSEA 富集曲线函数
gsea_curve <- function(fc, geneset) {
  hits <- fc$gene %in% geneset
  n_hits <- sum(hits)
  n_miss <- nrow(fc) - n_hits
  nr <- sum(abs(fc$logFC[hits]))
  es <- cumsum(ifelse(hits, abs(fc$logFC)/nr, -1/n_miss))
  return(data.frame(rank = 1:nrow(fc), es = es))
}

# 画 4 个基因集的 GSEA 富集曲线
plot_list <- list()
for (sn in set_names) {
  gs <- sets$gene[sets$set == sn]
  gs <- gs[gs %in% fc$gene]
  if (length(gs) < 3) next
  curve <- gsea_curve(fc, gs)
  plot_list[[sn]] <- curve
}

# 合成一张图（4个子图）
library(gridExtra)
plots <- lapply(names(plot_list), function(sn) {
  ggplot(plot_list[[sn]], aes(x = rank, y = es)) +
    geom_line(color = "#c0504d", linewidth = 0.8) +
    geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") +
    labs(title = sn, x = "Rank (high → low logFC)", y = "Enrichment Score") +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(size = 9, face = "bold"))
})

pdf(file.path(outdir, "planA_gsea_curves.pdf"), width = 9, height = 7)
grid.arrange(grobs = plots, ncol = 2)
dev.off()

# 单张核心图：PDX1 靶基因 GSEA 富集曲线（含命中基因标记）
gs <- sets$gene[sets$set == "PDX1_targets"]
gs <- gs[gs %in% fc$gene]
curve <- gsea_curve(fc, gs)
# 找 ES 最低点
min_es <- min(curve$es)
min_rank <- which.min(curve$es)

pdf(file.path(outdir, "planA_gsea_pdx1.pdf"), width = 7, height = 4.5)
p <- ggplot(curve, aes(x = rank, y = es)) +
  geom_line(color = "#c0504d", linewidth = 1) +
  geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") +
  geom_vline(xintercept = min_rank, color = "grey50", linetype = "dotted") +
  annotate("text", x = min_rank, y = min_es, label = sprintf("ES = %.3f", min_es),
           hjust = 0, vjust = 1.5, color = "#c0504d", fontface = "bold") +
  labs(title = "GSEA: PDX1 target genes enriched in downregulation (T2D)",
       subtitle = "PDX1 靶基因显著富集在 logFC 负端（下调）",
       x = "Rank (high → low logFC)", y = "Enrichment Score") +
  theme_bw(base_size = 11)
print(p)
dev.off()

cat("done\n")
