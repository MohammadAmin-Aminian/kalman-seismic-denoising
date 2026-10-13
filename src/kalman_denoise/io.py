"""Configuration and data loading."""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
from scipy.io import loadmat

__all__ = ["repo_root", "load_config", "load_noise", "load_earthquake"]


def repo_root():
    """Repository root (override with the KALMAN_DENOISE_ROOT environment variable)."""
    env = os.environ.get("KALMAN_DENOISE_ROOT")
    if env:
        return Path(env)
    here = Path(__file__).resolve()
    for p in here.parents:
        if (p / "config" / "default.json").is_file():
            return p
    return Path.cwd()


def load_config(path=None):
    path = Path(path) if path else repo_root() / "config" / "default.json"
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def _resolve(rel):
    p = Path(rel)
    return p if p.is_absolute() else repo_root() / p


def load_noise(rel_path):
    """Unit-variance Gaussian noise field stored as single precision (identical in MATLAB)."""
    return loadmat(_resolve(rel_path))["noise"].astype(float)


def load_earthquake(rel_path):
    m = loadmat(_resolve(rel_path), squeeze_me=True)
    return {
        "counts": np.asarray(m["counts"], dtype=float),
        "fs_hz": float(m["fs_hz"]),
        "station": str(m["station"]),
        "start_time": str(m["start_time"]),
    }
