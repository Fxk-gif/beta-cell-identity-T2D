"""Shared path resolution for the figure scripts.

The analysis data (r3_inputs/, sc_data/, ml_output/) live in the analysis
workspace, not inside this figures folder. Each script resolves the workspace
via the T2D_WORKSPACE environment variable when set; otherwise it walks upward
from this file's directory to find the folder containing r3_inputs/.
"""

import os


def find_workspace(start_dir=None):
    env = os.environ.get("T2D_WORKSPACE", "").rstrip("/\\")
    if env and os.path.isdir(env):
        return env

    d = os.path.abspath(start_dir or os.path.dirname(os.path.abspath(__file__)))
    for _ in range(12):
        if os.path.isdir(os.path.join(d, "r3_inputs")):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent

    raise RuntimeError(
        "Could not locate the analysis workspace (r3_inputs/ not found). "
        "Set the T2D_WORKSPACE environment variable."
    )
