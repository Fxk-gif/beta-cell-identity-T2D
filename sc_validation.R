# 方案A单细胞交叉验证：Xin et al. GSE81608
# 比较 PDX1/NKX6-1/SLC2A2/GLP1R 等β身份基因在 T2D vs ND β细胞中的表达
# 输入：sc_data/GSE81608_human_islets_rpkm.txt.gz（列=GSM） + GSE81608_annot.csv
library(ggplot2)
library(pROC)

sc_dir <- "C:/Users/45902/WorkBuddy/2026-09-20-15-12-44/sc_data"
outdir <- "C:/Users/45902/WorkBuddy/2026-09-20-15-12-44/figures"

# ---- 读注释 ----
annot <- read.csv(file.path(sc_dir, "GSE81608_annot.csv"), stringsAsFactors = FALSE)
annot$condition2 <- ifelse(annot$condition == "non-diabetic", "ND", "T2D")
beta <- annot[annot$subtype == "beta", ]
cat("beta 细胞数:", nrow(beta), "ND:", sum(beta$condition2=="ND"), "T2D:", sum(beta$condition2=="T2D"), "\n")

# ---- 读 rpkm 矩阵（gzip） ----
mat_file <- file.path(sc_dir, "GSE81608_human_islets_rpkm.txt.gz")
cat("读取 rpkm 矩阵...\n")
mat <- read.delim(gzfile(mat_file), check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)
cat("矩阵维度:", nrow(mat), "基因 x", ncol(mat), "细胞\n")

# 列名可能是 GSM 或去掉 GSM 的序号，检查
cols <- colnames(mat)
cat("列名示例:", paste(head(cols), collapse=", "), "\n")
# 匹配 beta 细胞的列
beta_cols <- intersect(beta$GSM, cols)
cat("匹配到的 beta 列数:", length(beta_cols), "（beta 注释总数", nrow(beta), "）\n")
if (length(beta_cols) < nrow(beta)) {
  # 尝试去掉 GSM 前缀匹配
  cols_short <- gsub("^GSM", "", cols)
  beta_short <- gsub("^GSM", "", beta$GSM)
  m <- match(beta_short, cols_short)
  cat("用去前缀方式匹配，命中:", sum(!is.na(m)), "\n")
  beta_cols <- cols[na.omit(m)]
}

# ---- 身份基因列表 ----
identity_genes <- c("PDX1","NKX6-1","SLC2A2","GLP1R","HNF1A","HNF4A","MAFA","MAFB",
                    "NEUROD1","PAX6","GCK","ABCC8","KCNJ11","IAPP","SLC30A8","INS")
neuro_genes <- c("GAD1","GAD2","CALCA","TSHR","NHLH2")
# 只保留矩阵中有的
identity_genes <- intersect(identity_genes, rownames(mat))
neuro_genes <- intersect(neuro_genes, rownames(mat))
cat("矩阵中的身份基因:", paste(identity_genes, collapse=", "), "\n")
cat("矩阵中的神经内分泌基因:", paste(neuro_genes, collapse=", "), "\n")

# ---- 提取 beta 细胞表达 ----
beta_expr <- mat[, beta_cols, drop = FALSE]
# 对齐分组
beta_order <- match(beta_cols, beta$GSM)
if (any(is.na(beta_order))) {
  beta_order <- match(gsub("^GSM","",beta_cols), gsub("^GSM","",beta$GSM))
}
beta_meta <- beta[beta_order, ]
stopifnot(nrow(beta_meta) == ncol(beta_expr))
beta_meta$condition2 <- factor(beta_meta$condition2, levels = c("ND","T2D"))

# ---- 身份评分（beta细胞层面） ----
zscore <- function(x) (x - mean(x)) / sd(x)
id_z <- t(apply(beta_expr[identity_genes, , drop=FALSE], 1, zscore))
identity_score <- colMeans(id_z, na.rm=TRUE)
beta_meta$identity_score <- identity_score

wi <- wilcox.test(identity_score ~ condition2, data = beta_meta)
cat("\n=== β细胞身份评分（单细胞） ===\n")
cat("Wilcoxon p =", wi$p.value, "\n")
cat("ND 均值:", round(mean(beta_meta$identity_score[beta_meta$condition2=="ND"]),3),
    "T2D 均值:", round(mean(beta_meta$identity_score[beta_meta$condition2=="T2D"]),3), "\n")

# ---- 每个身份基因的 T2D vs ND 比较 ----
gene_test <- data.frame(gene=character(), p=numeric(), nd_mean=numeric(), t2d_mean=numeric(),
                        log2fc=numeric(), stringsAsFactors=FALSE)
for (g in c(identity_genes, neuro_genes)) {
  vals <- as.numeric(beta_expr[g, ])
  nd <- vals[beta_meta$condition2=="ND"]
  t2d <- vals[beta_meta$condition2=="T2D"]
  p <- wilcox.test(nd, t2d)$p.value
  nd_m <- mean(nd); t2d_m <- mean(t2d)
  # log2fc（用 log1p 避免0）
  lfc <- log2((t2d_m + 0.1) / (nd_m + 0.1))
  gene_test <- rbind(gene_test, data.frame(gene=g, p=p, nd_mean=nd_m, t2d_mean=t2d_m, log2fc=lfc))
}
gene_test$padj <- p.adjust(gene_test$p, method="BH")
gene_test <- gene_test[order(gene_test$p), ]
cat("\n=== 单基因 T2D vs ND（β细胞）===\n")
print(gene_test, row.names=FALSE)

# ---- 图1：身份评分箱式图（单细胞） ----
p1 <- ggplot(beta_meta, aes(x=condition2, y=identity_score, fill=condition2)) +
  geom_boxplot(outlier.shape=NA, alpha=0.7) +
  geom_jitter(width=0.2, size=0.6, alpha=0.4) +
  scale_fill_manual(values=c("ND"="#4DBBD5","T2D"="#E64B35")) +
  labs(x=NULL, y="Beta cell identity score (z-scored)",
       title=sprintf("scRNA-seq (Xin et al.): β-cell identity loss in T2D (p=%.2e)", wi$p.value)) +
  theme_bw(base_size=12) + theme(legend.position="none")
ggsave(file.path(outdir, "planA_sc_identity_score.pdf"), p1, width=5.5, height=5)

# ---- 图2：关键身份基因点图（均值+显著性） ----
plot_genes <- intersect(c("PDX1","NKX6-1","SLC2A2","GLP1R","HNF1A","MAFA","NEUROD1","GCK","INS",
                          "GAD1","NHLH2","TSHR"), rownames(beta_expr))
plot_df <- data.frame()
for (g in plot_genes) {
  vals <- as.numeric(beta_expr[g, ])
  nd <- vals[beta_meta$condition2=="ND"]
  t2d <- vals[beta_meta$condition2=="T2D"]
  plot_df <- rbind(plot_df, data.frame(gene=g, ND=mean(nd), T2D=mean(t2d),
                                       p=wilcox.test(nd,t2d)$p.value))
}
library(tidyr)
plot_df_long <- pivot_longer(plot_df, cols=c(ND,T2D), names_to="group", values_to="expr")
plot_df_long$group <- factor(plot_df_long$group, levels=c("ND","T2D"))
plot_df_long$gene <- factor(plot_df_long$gene, levels=plot_genes)
plot_df$sig <- ifelse(plot_df$p<0.05, "*", ifelse(plot_df$p<0.01, "**", ""))
plot_df$sig <- ifelse(plot_df$p<0.001, "***", plot_df$sig)

p2 <- ggplot(plot_df_long, aes(x=gene, y=expr, fill=group)) +
  geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.6) +
  geom_errorbar(aes(ymin=expr, ymax=expr), width=0, position=position_dodge(0.7)) +
  scale_fill_manual(values=c("ND"="#4DBBD5","T2D"="#E64B35")) +
  geom_text(data=plot_df, aes(x=gene, y=max(plot_df_long$expr)*1.05, label=sig), inherit.aes=FALSE) +
  labs(x=NULL, y="Mean RPKM", fill=NULL, title="Identity & neuroendocrine genes in T2D vs ND β-cells (scRNA-seq)") +
  theme_bw(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(outdir, "planA_sc_genes.pdf"), p2, width=9, height=5)

# ---- 图3：单细胞身份评分 ROC（T2D vs ND） ----
outcome <- ifelse(beta_meta$condition2=="T2D", 1, 0)
roc_id <- roc(outcome, -identity_score, quiet=TRUE)
cat("\n=== 身份评分 ROC（单细胞）===\n")
cat("AUC =", round(as.numeric(roc_id$auc), 3), "\n")

# 单基因 ROC
roc_genes <- intersect(c("SLC2A2","PDX1","NKX6-1","GLP1R","HNF1A","MAFA"), rownames(beta_expr))
roc_list <- list("Identity score"=roc_id)
for (g in roc_genes) {
  vals <- as.numeric(beta_expr[g, ])
  roc_list[[g]] <- roc(outcome, -vals, quiet=TRUE)
}
aucs <- sapply(roc_list, function(x) as.numeric(x$auc))
pdf(file.path(outdir, "planA_sc_roc.pdf"), width=6.5, height=6.5)
cols <- rainbow(length(roc_list))
plot(roc_list[[1]], col=cols[1], lwd=2, main="ROC (scRNA-seq β-cells, Xin et al.)")
for (i in 2:length(roc_list)) plot(roc_list[[i]], col=cols[i], lwd=2, add=TRUE)
legend("bottomright", legend=sprintf("%s (AUC=%.3f)", names(aucs), aucs), col=cols, lwd=2, cex=0.7, bty="n")
dev.off()
auc_df <- data.frame(gene=names(aucs), AUC=round(aucs,3))
write.csv(auc_df, file.path(outdir, "planA_sc_roc_auc.csv"), row.names=FALSE)
cat("\n=== 单细胞 ROC AUC ===\n")
print(auc_df[order(-auc_df$AUC),], row.names=FALSE)

cat("\n=== 单细胞交叉验证完成 ===\n")
