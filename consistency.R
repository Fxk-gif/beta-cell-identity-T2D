# 方案A 核心证据图：β身份评分跨数据集一致性
library(ggplot2)
library(tidyr)

base <- "C:/Users/45902/WorkBuddy/2026-09-20-15-12-44/r3_inputs"
outdir <- "C:/Users/45902/WorkBuddy/2026-09-20-15-12-44/figures"

zscore <- function(x) (x - mean(x)) / sd(x)

calc_score <- function(ds) {
  expr <- read.csv(file.path(base, paste0(ds, "_expr.csv")), check.names=FALSE, row.names=1)
  grp  <- read.csv(file.path(base, paste0(ds, "_groups.csv")), stringsAsFactors=FALSE)
  colnames(grp) <- c("sample","group")
  grp <- grp[grp$sample %in% colnames(expr), ]
  expr <- expr[, grp$sample]
  keep <- grp$group %in% c("T2D","ND")
  expr <- expr[, keep]; grp <- grp[keep,]
  beta_identity <- intersect(c("PDX1","NKX6-1","HNF1A","HNF4A","MAFB","NEUROD1","PAX6",
        "SLC2A2","GCK","ABCC8","IAPP","GLP1R","SLC30A8","PCSK1","PCSK2","INS"), rownames(expr))
  beta_z <- t(apply(expr[beta_identity,,drop=FALSE],1,zscore))
  data.frame(dataset=ds, group=grp$group,
             identity_score=colMeans(beta_z,na.rm=TRUE))
}

tr <- calc_score("GSE76895")
va <- calc_score("GSE76894")
dat <- rbind(tr, va)
dat$dataset <- factor(dat$dataset, levels=c("GSE76895","GSE76894"),
                      labels=c("Discovery (GSE76895)","Validation (GSE76894)"))
dat$group <- factor(dat$group, levels=c("ND","T2D"))

p <- ggplot(dat, aes(x=group, y=identity_score, fill=group)) +
  geom_boxplot(outlier.shape=NA, alpha=0.7) +
  geom_jitter(width=0.2, size=1.2, alpha=0.5) +
  facet_wrap(~dataset, scales="free_y") +
  scale_fill_manual(values=c("ND"="#4DBBD5","T2D"="#E64B35")) +
  labs(x=NULL, y="Beta cell identity score (z-scored)",
       title="Beta cell identity loss is reproducible across two independent datasets") +
  theme_bw(base_size=12) + theme(legend.position="none",
        plot.title=element_text(size=11, hjust=0.5))
ggsave(file.path(outdir, "planA_identity_consistency.pdf"), p, width=8, height=4.2)

cat("跨数据集一致性图已输出\n")
