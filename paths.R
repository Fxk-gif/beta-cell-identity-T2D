# unified path configuration (reproducible) -- R version
# usage: run R from the repository root, then source("paths.R")
# data-directory order: env var T2D_REPO_DATA, else data/ under the working directory

data_root <- Sys.getenv("T2D_REPO_DATA", unset = NA)
if (is.na(data_root) || nchar(data_root) == 0) {
  data_root <- file.path(getwd(), "data")
}

r3_dir   <- file.path(data_root, "r3_inputs")
sc_dir   <- file.path(data_root, "sc_data")
fig_dir  <- file.path(data_root, "figures")
geo_proc <- file.path(data_root, "geo_data", "processed")

for (d in c(r3_dir, sc_dir, fig_dir, geo_proc)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

cat("data root:", data_root, "\n")
