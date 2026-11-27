# -*- coding: utf-8 -*-
"""
Unified path configuration (reproducible).

All analysis scripts locate the data directories through this module, avoiding hard-coded absolute local paths.

Data-directory resolution order:
  1. environment variable `T2D_REPO_DATA` (if set)
  2. the `data/` subdirectory under the repository root

Usage:
    import paths
    df = pd.read_csv(paths.r3_inputs('GSE76895_expr.csv'), index_col=0)

Directory layout (place the raw GEO data in the corresponding subdirectories, or run build_r3_inputs.py to generate them):
    data/
      geo_data/processed/   # processed bulk microarray expression matrices and groups
      r3_inputs/            # analysis inputs (GSE76895/GSE76894_expr.csv, *_groups.csv)
      sc_data/              # single-cell expression and annotations (xin_beta_expr.csv etc.)
      ml_output/            # model outputs
      figures/              # figure outputs
"""
import os

_ROOT = os.path.dirname(os.path.abspath(__file__))

# data root: prefer the environment variable, else data/ under the repo
_DATA = os.environ.get('T2D_REPO_DATA', os.path.join(_ROOT, 'data'))


def _join(*parts):
    return os.path.join(_DATA, *parts)


def geo_processed(*parts):
    """Directory of the processed bulk microarray data."""
    return _join('geo_data', 'processed', *parts)


def r3_inputs(*parts):
    """Analysis-input directory."""
    return _join('r3_inputs', *parts)


def sc_data(*parts):
    """Single-cell data directory."""
    return _join('sc_data', *parts)


def ml_output(*parts):
    """Model-output directory."""
    return _join('ml_output', *parts)


def figures(*parts):
    """Figure-output directory."""
    return _join('figures', *parts)


def ensure_dirs():
    """Create an output directory."""
    for d in (ml_output(), figures(), figures('png')):
        os.makedirs(d, exist_ok=True)


if __name__ == '__main__':
    print('Data root:', _DATA)
    print('r3_inputs:', r3_inputs())
    print('sc_data:', sc_data())
