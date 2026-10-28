# 方案A验证集交叉验证：GSE76894 独立验证
library(pROC)
library(ggplot2)
library(tidyr)

base <- "C:/Users/45902/WorkBuddy/2026-09-20-15-12-44/r3_inputs"
outdir <- "C:/Users/45902/WorkBuddy/2026-09-20-15-12-44/figures"

expr <- read.csv(file.path(base, "GSE76894_expr.csv"), check.names = FALSE, row.names = 1)
grp  <- read.csv(file.path(base, "GSE76894_groups.csv"), stringsAsFactors = FALSE)
colnames(grp) <- c("sample", "group")
grp <- grp[grp$sample %in% colnames(expr), ]
expr <- expr[, grp$sample]
keep <- grp$group %in% c("T2D", "ND")
expr <- expr[, keep]
grp  <- grp[keep, ]
cat("验证集维度:", nrow(expr), "基因 x", ncol(expr), "样本\n")
cat("分组:", table(grp$group), "\n")

beta_identity <- c("PDX1","NKX6-1","HNF1A","HNF4A","MAFB","NEUROD1","PAX6",
                   "SLC2A2","GCK","ABCC8","IAPP","GLP1R","SLC30A8","PCSK1","PCSK2","INS")
neuro_trans <- c("GAD1","GAD2","CALCA","TSHR","NHLH2")
beta_identity <- intersect(beta_identity, rownames(expr))
neuro_trans   <- intersect(neuro_trans, rownames(expr))

zscore <- function(x) (x - mean(x)) / sd(x)
beta_z   <- t(apply(expr[beta_identity, , drop=FALSE], 1, zscore))
identity_score <- colMeans(beta_z, na.rm=TRUE)
neuro_z  <- t(apply(expr[neuro_trans, , drop=FALSE], 1, zscore))
neuro_score <- colMeans(neuro_z, na.rm=TRUE)

df <- data.frame(group = grp$group, identity_score, neuro_score)
df$group <- factor(df$group, levels=c("ND","T2D"))
wi <- wilcox.test(identity_score ~ group, data=df)
wn <- wilcox.test(neuro_score ~ group, data=df)
cat("\n=== 验证集身份评分 ===\n")
cat("β身份评分 p =", wi$p.value, "\n")
cat("神经内分泌评分 p =", wn$p.value, "\n")
cat("β身份评分均值 ND vs T2D:", round(mean(df$identity_score[df$group=="ND"]),3),
    "vs", round(mean(df$identity_score[df$group=="T2D"]),3), "\n")
cat("神经内分泌评分均值 ND vs T2D:", round(mean(df$neuro_score[df$group=="ND"]),3),
    "vs", round(mean(df$neuro_score[df$group=="T2D"]),3), "\n")

# 关键基因ROC（验证集）
roc_genes <- intersect(c("SLC2A2","PDX1","NKX6-1","HNF1A","GLP1R"), rownames(expr))
neuro_up  <- intersect(c("GAD1","NHLH2","TSHR","CALCA"), rownames(expr))
outcome <- ifelse(df$group=="T2D", 1, 0)
roc_list <- list()
for (g in roc_genes) roc_list[[paste0(g," (down)")]] <- roc(outcome, -as.numeric(expr[g,]), quiet=TRUE)
for (g in neuro_up)  roc_list[[paste0(g," (up)")]]   <- roc(outcome,  as.numeric(expr[g,]), quiet=TRUE)

aucs <- sapply(roc_list, function(x) as.numeric(x$auc))
auc_df <- data.frame(gene=names(aucs), AUC=round(aucs,3))
write.csv(auc_df, file.path(outdir, "planA_roc_auc_val.csv"), row.names=FALSE)
cat("\n=== 验证集身份基因 ROC AUC ===\n")
print(auc_df[order(-auc_df$AUC),])

pdf(file.path(outdir, "planA_roc_val.pdf"), width=7, height=7)
cols <- rainbow(length(roc_list))
plot(roc_list[[1]], col=cols[1], lwd=2, main="ROC (validation, GSE76894)")
for (i in 2:length(roc_list)) plot(roc_list[[i]], col=cols[i], lwd=2, add=TRUE)
legend("bottomright", legend=sprintf("%s (AUC=%.3f)", names(aucs), aucs), col=cols, lwd=2, cex=0.7, bty="n")
dev.off()

# 评分箱式图（验证集）
df_long <- pivot_longer(df, cols=c(identity_score, neuro_score), names_to="score", values_to="value")
df_long$score <- factor(df_long$score, levels=c("identity_score","neuro_score"),
                        labels=c("Beta cell identity score","Neuroendocrine score"))
p <- ggplot(df_long, aes(x=group, y=value, fill=group)) +
  geom_boxplot(outlier.shape=NA, alpha=0.7) + geom_jitter(width=0.2, size=1.2, alpha=0.5) +
  facet_wrap(~score, scales="free_y") +
  scale_fill_manual(values=c("ND"="#4DBBD5","T2D"="#E64B35")) +
  labs(x=NULL, y="Standardized score", title="Validation (GSE76894): identity loss & neuroendocrine shift") +
  theme_bw(base_size=12) + theme(legend.position="none")
ggsave(file.path(outdir, "planA_score_box_val.pdf"), p, width=7, height=4)

cat("\n=== 验证集分析完成 ===\n")
