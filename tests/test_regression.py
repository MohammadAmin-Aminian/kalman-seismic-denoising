"""The published numbers (results/metrics.json) must be reproduced exactly."""
import json

import numpy as np

from kalman_denoise import pipeline as P
from kalman_denoise.io import load_config, repo_root


def _published():
    with open(repo_root() / "results" / "metrics.json", encoding="utf-8") as fh:
        return json.load(fh)


def test_synthetic_snr_table():
    cfg, pub = load_config(), _published()
    for s in cfg["noise"]["sigma_levels"]:
        r = P.run_synthetic(cfg, s)
        for m, v in pub["synthetic_snr_db"][f"{s:g}"].items():
            assert abs(r["snr"][m] - v) < 1e-6, (s, m)


def test_chirp_matches_theory():
    ch = P.run_chirp(load_config())
    assert ch["max_gain_error"] < 0.05


def test_earthquake_filter_numbers():
    eq, pub = P.run_earthquake(load_config()), _published()["earthquake"]
    for k in ("main", "alt"):
        assert np.isclose(eq[k]["K"], pub[k]["K"]) and np.isclose(eq[k]["rms_ratio"], pub[k]["rms_ratio"])
