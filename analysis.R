# 方案A：β细胞身份丢失/神经内分泌转分化 核心分析
# 输出到 figures/planA_* 

library(pROC)
library(ggplot2)

base <- "C:/Users/45902/WorkBuddy/2026-09-20-15-12-44/r3_inputs"
outdir <- "C:/Users/45902/WorkBuddy/2026-09-20-15-12-44/figures"
dir.create(outdir, showWarnings = FALSE)

# ---- 读取训练集表达矩阵 ----
expr <- read.csv(file.path(base, "GSE76895_expr.csv"), check.names = FALSE, row.names = 1)
grp  <- read.csv(file.path(base, "GSE76895_groups.csv"), stringsAsFactors = FALSE)
colnames(grp) <- c("sample", "group")
grp <- grp[grp$sample %in% colnames(expr), ]
expr <- expr[, grp$sample]
# 只保留 T2D vs ND
keep <- grp$group %in% c("T2D", "ND")
expr <- expr[, keep]
grp  <- grp[keep, ]
cat("训练集维度:", nrow(expr), "基因 x", ncol(expr), "样本\n")
cat("分组:", table(grp$group), "\n")

# ---- 基因集定义 ----
beta_identity <- c("PDX1","NKX6-1","HNF1A","HNF4A","MAFA","MAFB","NEUROD1","PAX6",
                   "SLC2A2","GCK","ABCC8","KCNJ11","IAPP","UCN3","GLP1R","SLC30A8","PCSK1","PCSK2","INS")
neuro_trans <- c("GAD1","GAD2","CALCA","TSHR","NHLH2")
dediff     <- c("FOXE1","MYCN","CHL1","OPRD1","GRIN2D","ALDH1A3","HES1","SOX9","NGN3")

# 只保留有探针的基因
beta_identity <- intersect(beta_identity, rownames(expr))
neuro_trans   <- intersect(neuro_trans, rownames(expr))
dediff        <- intersect(dediff, rownames(expr))
cat("有探针的 beta_identity:", paste(beta_identity, collapse=", "), "\n")
cat("有探针的 neuro_trans:", paste(neuro_trans, collapse=", "), "\n")
cat("有探针的 dediff:", paste(dediff, collapse=", "), "\n")

# ---- β细胞身份评分（z-score 平均） ----
zscore <- function(x) (x - mean(x)) / sd(x)
beta_mat <- expr[beta_identity, , drop = FALSE]
beta_z   <- t(apply(beta_mat, 1, zscore))
identity_score <- colMeans(beta_z, na.rm = TRUE)

# 神经内分泌评分
neuro_mat <- expr[neuro_trans, , drop = FALSE]
neuro_z   <- t(apply(neuro_mat, 1, zscore))
neuro_score <- colMeans(neuro_z, na.rm = TRUE)

# ---- 身份评分组间比较 ----
df <- data.frame(
  sample = grp$sample,
  group = grp$group,
  identity_score = identity_score,
  neuro_score = neuro_score,
  stringsAsFactors = FALSE
)
df$group <- factor(df$group, levels = c("ND", "T2D"))

# 统计检验
wi <- wilcox.test(identity_score ~ group, data = df)
wn <- wilcox.test(neuro_score ~ group, data = df)
cat("\n=== 身份评分 Wilcoxon 检验 ===\n")
cat("β身份评分: p =", wi$p.value, "\n")
cat("神经内分泌评分: p =", wn$p.value, "\n")
cat("β身份评分均值 ND vs T2D:", 
    round(mean(df$identity_score[df$group=="ND"]),3), "vs",
    round(mean(df$identity_score[df$group=="T2D"]),3), "\n")
cat("神经内分泌评分均值 ND vs T2D:",
    round(mean(df$neuro_score[df$group=="ND"]),3), "vs",
    round(mean(df$neuro_score[df$group=="T2D"]),3), "\n")

# ---- 图1：身份评分 + 神经内分泌评分 箱式图 ----
library(tidyr)
df_long <- pivot_longer(df, cols = c(identity_score, neuro_score),
                        names_to = "score", values_to = "value")
df_long$score <- factor(df_long$score, 
                        levels = c("identity_score", "neuro_score"),
                        labels = c("β cell identity score", "Neuroendocrine score"))

p1 <- ggplot(df_long, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.5) +
  facet_wrap(~ score, scales = "free_y") +
  scale_fill_manual(values = c("ND" = "#4DBBD5", "T2D" = "#E64B35")) +
  labs(x = NULL, y = "Standardized score", title = "β cell identity loss & neuroendocrine trans-differentiation in T2D islets") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
ggsave(file.path(outdir, "planA_score_box.pdf"), p1, width = 7, height = 4)

# ---- 图2：关键身份基因表达热图（T2D vs ND 各基因均值） ----
key_genes <- c(beta_identity, neuro_trans, dediff)
key_genes <- unique(key_genes)
mean_expr <- t(apply(expr[key_genes, , drop = FALSE], 1, function(r) tapply(r, grp$group, mean)))
mean_expr <- data.frame(mean_expr)
mean_expr$gene <- rownames(mean_expr)
mean_expr$category <- ifelse(mean_expr$gene %in% beta_identity, "β identity",
                      ifelse(mean_expr$gene %in% neuro_trans, "Neuroendocrine", "Dediff/Progenitor"))
mean_long <- pivot_longer(mean_expr, cols = c("ND", "T2D"), names_to = "group", values_to = "expr")
mean_long$group <- factor(mean_long$group, levels = c("ND", "T2D"))
mean_long$gene <- factor(mean_long$gene, levels = rev(key_genes))

p2 <- ggplot(mean_long, aes(x = group, y = gene, fill = expr)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#E64B35", midpoint = median(mean_long$expr)) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL, fill = "Mean expr", title = "Mean expression of identity-related genes (T2D vs ND)") +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_text(size = 8))
ggsave(file.path(outdir, "planA_identity_heatmap.pdf"), p2, width = 6, height = 8)

# ---- 图3：关键基因 ROC（身份丢失基因用反向预测 T2D） ----
# 身份基因在T2D下调 → 用 (1 - zscore) 或直接反向 做ROC
roc_genes <- c("SLC2A2","PDX1","NKX6-1","HNF1A","GLP1R")  # β身份下调基因
roc_genes <- intersect(roc_genes, rownames(expr))
neuro_up <- c("GAD1","NHLH2","TSHR","CALCA")  # 神经内分泌上调基因
neuro_up <- intersect(neuro_up, rownames(expr))

outcome <- ifelse(df$group == "T2D", 1, 0)
roc_list <- list()
for (g in roc_genes) {
  vals <- as.numeric(expr[g, ])
  roc_list[[paste0(g, " (down)")]] <- roc(outcome, -vals, quiet = TRUE)  # 下调→反向
}
for (g in neuro_up) {
  vals <- as.numeric(expr[g, ])
  roc_list[[paste0(g, " (up)")]] <- roc(outcome, vals, quiet = TRUE)
}

pdf(file.path(outdir, "planA_roc.pdf"), width = 7, height = 7)
cols <- rainbow(length(roc_list))
aucs <- sapply(roc_list, function(x) as.numeric(x$auc))
names(aucs) <- names(roc_list)
plot(roc_list[[1]], col = cols[1], lwd = 2, main = "ROC - identity-loss & neuroendocrine genes")
for (i in 2:length(roc_list)) {
  plot(roc_list[[i]], col = cols[i], lwd = 2, add = TRUE)
}
legend("bottomright", legend = sprintf("%s (AUC=%.3f)", names(aucs), aucs),
       col = cols, lwd = 2, cex = 0.7, bty = "n")
dev.off()

# 输出AUC表
auc_df <- data.frame(gene = names(aucs), AUC = round(aucs, 3))
write.csv(auc_df, file.path(outdir, "planA_roc_auc.csv"), row.names = FALSE)
cat("\n=== 身份基因 ROC AUC ===\n")
print(auc_df[order(-auc_df$AUC), ])

cat("\n=== 方案A核心分析完成，图已输出到 figures/planA_* ===\n")
