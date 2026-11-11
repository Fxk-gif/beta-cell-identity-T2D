# Figure 3. Diagnostic performance of individual beta-cell identity genes
#
# Panels:
#   (A) ROC curves for individual identity genes in the discovery dataset (GSE76895)
#   (B) ROC curves in the validation dataset (GSE76894)
#
# The five identity genes most consistently downregulated in T2DM islets are
# tested for their ability to discriminate T2D from ND samples. Because these
# genes are reduced in T2D, the predictor is the negated expression so that a
# higher value indicates T2D. AUC is computed with the pROC package.

suppressPackageStartupMessages({
  library(pROC)
  library(ggplot2)
  library(gridExtra)
})

# ---- resolve the analysis workspace -----------------------------------------
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

# ---- gene selection ---------------------------------------------------------
# Five identity genes reduced in T2D, used as discriminators (negated).
IDENTITY_DOWN <- c("SLC2A2", "PDX1", "NKX6-1", "HNF1A", "GLP1R")

ROC_COLORS <- c("SLC2A2" = "#E64B35", "PDX1" = "#3B6FB6",
                "NKX6-1" = "#4DBBD5", "HNF1A" = "#2F9E6F",
                "GLP1R" = "#C792EA")

# ---- load data and compute ROC ----------------------------------------------
build_roc <- function(ds) {
  expr <- read.csv(file.path(DATA, paste0(ds, "_expr.csv")),
                   check.names = FALSE, row.names = 1)
  grp  <- read.csv(file.path(DATA, paste0(ds, "_groups.csv")), header = FALSE)
  colnames(grp) <- c("sample", "group")
  grp  <- grp[grp$group %in% c("T2D", "ND"), ]
  expr <- expr[, grp$sample, drop = FALSE]

  outcome <- ifelse(grp$group == "T2D", 1, 0)
  curves <- list()
  aucs   <- numeric()
  for (g in IDENTITY_DOWN) {
    if (!(g %in% rownames(expr))) next
    vals <- as.numeric(expr[g, ])
    r <- roc(outcome, -vals, quiet = TRUE, direction = "<")
    curves[[g]] <- r
    aucs[g] <- as.numeric(r$auc)
  }
  list(curves = curves, aucs = aucs)
}

r95 <- build_roc("GSE76895")
r94 <- build_roc("GSE76894")

cat("Discovery AUC:\n");  print(round(r95$aucs, 3))
cat("Validation AUC:\n"); print(round(r94$aucs, 3))

# ---- plot ---------------------------------------------------------------
roc_plot <- function(roc_list, aucs, title) {
  plot_data <- do.call(rbind, lapply(names(roc_list), function(g) {
    r <- roc_list[[g]]
    data.frame(gene = g, fpr = as.numeric(1 - r$specificities),
               tpr = as.numeric(r$sensitivities), auc = aucs[g],
               row.names = NULL)
  }))
  # keep curve order stable (as listed)
  plot_data$gene <- factor(plot_data$gene, levels = names(roc_list))

  ggplot(plot_data, aes(x = fpr, y = tpr, color = gene)) +
    geom_abline(intercept = 0, slope = 1, color = "grey70", linetype = "dashed") +
    geom_line(linewidth = 1) +
    scale_color_manual(
      values = ROC_COLORS[names(roc_list)],
      labels = sprintf("%s (AUC = %.3f)", names(roc_list), aucs[names(roc_list)])
    ) +
    labs(x = "1 - Specificity", y = "Sensitivity", color = NULL) +
    ggtitle(title) +
    coord_equal() +
    theme_bw(base_size = 11) +
    theme(legend.position = c(0.75, 0.35),
          legend.text = element_text(size = 8),
          legend.background = element_rect(fill = "white", color = "grey80"),
          plot.title = element_text(size = 11, face = "bold"))
}

pA <- roc_plot(r95$curves, r95$aucs, "ROC - identity genes (discovery, GSE76895)")
pB <- roc_plot(r94$curves, r94$aucs, "ROC - identity genes (validation, GSE76894)")

# Figure 3 (single file, two panels side by side)
p3 <- grid.arrange(pA, pB, ncol = 2)
ggsave(file.path(OUT, "Figure3.TIF"), p3, width = 9.2, height = 4.6,
       dpi = 300, device = "tiff", compression = "lzw")

cat("Figure 3 written to", OUT, "\n")
