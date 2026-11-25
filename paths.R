# 统一路径配置（可复现）—— R 版本
# 用法：在仓库根目录下运行 R，然后 source("paths.R")
# 数据目录优先级：环境变量 T2D_REPO_DATA，否则当前工作目录下的 data/

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

cat("数据根目录:", data_root, "\n")
