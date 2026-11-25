# -*- coding: utf-8 -*-
"""
统一路径配置（可复现）。

所有分析脚本通过本模块定位数据目录，避免硬编码本机绝对路径。

数据目录优先级：
  1. 环境变量 `T2D_REPO_DATA`（若设置）
  2. 本仓库根目录下的 `data/` 子目录

用法：
    import paths
    df = pd.read_csv(paths.r3_inputs('GSE76895_expr.csv'), index_col=0)

目录结构约定（请先将 GEO 原始数据放到对应子目录，或运行 build_r3_inputs.py 生成）：
    data/
      geo_data/processed/   # bulk 微阵列处理后的表达矩阵与分组
      r3_inputs/            # 分析输入（GSE76895/GSE76894_expr.csv, *_groups.csv）
      sc_data/              # 单细胞表达与注释（xin_beta_expr.csv 等）
      ml_output/            # 模型输出
      figures/              # 图输出
"""
import os

_ROOT = os.path.dirname(os.path.abspath(__file__))

# 数据根目录：优先环境变量，否则仓库下 data/
_DATA = os.environ.get('T2D_REPO_DATA', os.path.join(_ROOT, 'data'))


def _join(*parts):
    return os.path.join(_DATA, *parts)


def geo_processed(*parts):
    """bulk 微阵列处理后的数据目录。"""
    return _join('geo_data', 'processed', *parts)


def r3_inputs(*parts):
    """分析输入目录。"""
    return _join('r3_inputs', *parts)


def sc_data(*parts):
    """单细胞数据目录。"""
    return _join('sc_data', *parts)


def ml_output(*parts):
    """模型输出目录。"""
    return _join('ml_output', *parts)


def figures(*parts):
    """图输出目录。"""
    return _join('figures', *parts)


def ensure_dirs():
    """创建输出目录。"""
    for d in (ml_output(), figures(), figures('png')):
        os.makedirs(d, exist_ok=True)


if __name__ == '__main__':
    print('数据根目录:', _DATA)
    print('r3_inputs:', r3_inputs())
    print('sc_data:', sc_data())
