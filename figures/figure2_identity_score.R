# Figure 2. Reduced beta-cell identity score in T2DM islets across two separate datasets
#
# Panels:
#   (A) discovery (GSE76895) identity-score box plot, T2D vs ND
#   (B) validation (GSE76894) identity-score box plot, T2D vs ND
#   (C) side-by-side comparison of the two cohorts
#
# The identity score is the per-sample mean of the z-scored expression of the
# 15-gene identity panel (the 16-gene panel minus PDX1, which has no probe on
# the GSE76894 platform). Group comparison uses the Wilcoxon rank-sum test with
# exact p-values (wilcox.test(..., exact = TRUE)).

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
})

# ---- resolve the analysis workspace (contains r3_inputs/) -------------------
# The data live in the analysis workspace, not in this script's folder. The
# location is taken from the T2D_WORKSPACE environment variable when set;
# otherwise the script walks upward from its own directory to find a folder
# that contains r3_inputs/.
find_workspace <- function(start) {
  env <- Sys.getenv("T2D_WORKSPACE", unset = "")
  if (nzchar(env) && dir.exists(env)) {
    return(env)
  }
  d <- normalizePath(start, mustWork = FALSE)
  for (i in seq_len(12)) {
    if (dir.exists(file.path(d, "r3_inputs"))) {
      return(d)
    }
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Could not locate the analysis workspace (r3_inputs/ not found).",
       "Set the T2D_WORKSPACE environment variable.")
}

WS   <- find_workspace(getwd())
DATA <- file.path(WS, "r3_inputs")
OUT  <- file.path(WS, "figures")
dir.create(OUT, showWarnings = FALSE)

# ---- gene panel -------------------------------------------------------------
TF    <- c("PDX1", "NKX6-1", "HNF1A", "HNF4A", "MAFB", "NEUROD1", "PAX6")
FUNC  <- c("SLC2A2", "GCK", "ABCC8", "IAPP", "GLP1R", "SLC30A8",
           "PCSK1", "PCSK2", "INS")
G16   <- c(TF, FUNC)          # 16-gene panel
G15   <- setdiff(G16, "PDX1") # 15-gene panel (PDX1 absent on GSE76894)

ND_COL  <- "#4DBBD5"
T2D_COL <- "#E64B35"

# ---- per-cohort identity score + exact Wilcoxon -----------------------------
compute_score <- function(ds) {
  expr <- read.csv(file.path(DATA, paste0(ds, "_expr.csv")),
                   check.names = FALSE, row.names = 1)
  grp  <- read.csv(file.path(DATA, paste0(ds, "_groups.csv")), header = FALSE)
  colnames(grp) <- c("sample", "group")
  grp  <- grp[grp$group %in% c("T2D", "ND"), ]
  expr <- expr[, grp$sample, drop = FALSE]

  genes <- intersect(G15, rownames(expr))
  sub   <- as.matrix(expr[genes, , drop = FALSE])
  storage.mode(sub) <- "double"
  z     <- t(apply(sub, 1, function(x) (x - mean(x, na.rm = TRUE)) /
                     sd(x, na.rm = TRUE)))
  score <- colMeans(z, na.rm = TRUE)

  nd  <- score[grp$group == "ND"]
  t2d <- score[grp$group == "T2D"]
  p   <- wilcox.test(t2d, nd, exact = TRUE)$p.value

  data.frame(sample = names(score), group = grp$group, score = score,
             n_genes = length(genes), p = p, stringsAsFactors = FALSE)
}

d95 <- compute_score("GSE76895")
d94 <- compute_score("GSE76894")
d95$cohort <- "Discovery (GSE76895)"
d94$cohort <- "Validation (GSE76894)"
dat <- rbind(d95, d94)
dat$group <- factor(dat$group, levels = c("ND", "T2D"))

# ---- plotting helper --------------------------------------------------------
p_label <- function(p) {
  if (p < 0.001) {
    return(sprintf("p = %.1e", p))
  }
  sprintf("p = %.3f", p)
}

box_plot <- function(df, cohort_name, show_label = FALSE) {
  sub <- df[df$cohort == cohort_name, ]
  pval <- sub$p[1]
  ymax <- max(sub$score)
  ymin <- min(sub$score)
  p <- ggplot(sub, aes(x = group, y = score, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.75) +
    geom_jitter(width = 0.12, size = 1.6, alpha = 0.75,
                color = "grey30") +
    scale_fill_manual(values = c("ND" = ND_COL, "T2D" = T2D_COL)) +
    labs(x = NULL, y = "Beta-cell identity score (z-scored)") +
    ggtitle(cohort_name) +
    ylim(ymin - 0.25 * (ymax - ymin), ymax + 0.25 * (ymax - ymin)) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(size = 11, face = "bold"))
  if (show_label) {
    p <- p + annotate("text", x = 1.5, y = ymax + 0.12 * (ymax - ymin),
                      label = p_label(pval), size = 3.5)
  }
  p
}

pA <- box_plot(dat, "Discovery (GSE76895)", show_label = TRUE)
pB <- box_plot(dat, "Validation (GSE76894)", show_label = TRUE)

cat(sprintf("Discovery: n_genes=%d  p(exact)=%.3e\n", d95$n_genes[1], d95$p[1]))
cat(sprintf("Validation: n_genes=%d  p(exact)=%.3e\n", d94$n_genes[1], d94$p[1]))

# Figure 2 (single file, two panels side by side)
p2 <- grid.arrange(pA, pB, ncol = 2,
                   top = "Reduced beta-cell identity score is reproducible across two separate datasets")
ggsave(file.path(OUT, "Figure2.TIF"), p2, width = 8.4, height = 3.8,
       dpi = 300, device = "tiff", compression = "lzw")

cat("Figure 2 written to", OUT, "\n")
