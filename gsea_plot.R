source("paths.R")
library(ggplot2)

base <- r3_dir
outdir <- fig_dir

# read the genome-wide logFC
fc <- read.csv(file.path(base, "volcano_input.csv"))
names(fc) <- c("gene", "logFC", "padj")
fc <- fc[order(-fc$logFC), ]  # descending

# read the target gene sets
sets <- read.csv(file.path(base, "tf_target_sets.csv"))
set_names <- unique(sets$set)

# GSEA enrichment-curve function
gsea_curve <- function(fc, geneset) {
  hits <- fc$gene %in% geneset
  n_hits <- sum(hits)
  n_miss <- nrow(fc) - n_hits
  nr <- sum(abs(fc$logFC[hits]))
  es <- cumsum(ifelse(hits, abs(fc$logFC)/nr, -1/n_miss))
  return(data.frame(rank = 1:nrow(fc), es = es))
}

# plot the GSEA enrichment curves of the 4 gene sets
plot_list <- list()
for (sn in set_names) {
  gs <- sets$gene[sets$set == sn]
  gs <- gs[gs %in% fc$gene]
  if (length(gs) < 3) next
  curve <- gsea_curve(fc, gs)
  plot_list[[sn]] <- curve
}

# combine into one figure (4 panels)
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

# single core panel: PDX1 target-gene GSEA enrichment curve (with hit markers)
gs <- sets$gene[sets$set == "PDX1_targets"]
gs <- gs[gs %in% fc$gene]
curve <- gsea_curve(fc, gs)
# find the ES minimum
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
       subtitle = "PDX1 targets significantly enriched at the negative (downregulated) logFC end",
       x = "Rank (high → low logFC)", y = "Enrichment Score") +
  theme_bw(base_size = 11)
print(p)
dev.off()

cat("done\n")
