"""Build the files in data/ (run once; the outputs are committed).

    python tools/make_data.py --sigs /path/to/original/Sigs.mat

* data/synthetic_noise.mat  unit-variance Gaussian noise (nt x ntraces), single precision
* data/chirp_noise.mat      unit-variance Gaussian noise for the chirp test, single precision
* data/earthquake_IU_DAV_BHZ.mat  vertical broadband record (counts) taken from Sigs.mat/sig1

The noise is drawn once with a fixed seed and stored, so that Python and
MATLAB use bit-identical inputs (their random generators differ).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from kalman_denoise.io import load_config  # noqa: E402
from kalman_denoise.synthetic import hyperbolic_gather, time_axis  # noqa: E402

SEED = 20210926  # date of the original synthetic experiment script


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sigs", help="original Sigs.mat (only needed for the earthquake file)")
    a = ap.parse_args()
    cfg = load_config()
    rng = np.random.default_rng(SEED)

    d, _, _ = hyperbolic_gather(cfg["synthetic"])
    noise = rng.standard_normal(d.shape).astype(np.float32)
    savemat(ROOT / cfg["noise"]["file"], {"noise": noise, "seed": SEED}, do_compression=True)

    c = cfg["chirp"]
    nt = time_axis(c["dt_s"], c["duration_s"]).size
    cn = rng.standard_normal((nt, 1)).astype(np.float32)
    savemat(ROOT / c["noise_file"], {"noise": cn, "seed": SEED}, do_compression=True)

    if a.sigs:
        s = loadmat(a.sigs)
        sig1 = np.asarray(s["sig1"]).ravel().astype(np.int32)
        savemat(
            ROOT / cfg["earthquake"]["file"],
            {
                "counts": sig1[:, None],
                "fs_hz": 40.0,
                "station": "IU.DAV.00.BHZ",
                "start_time": "2021-09-05T15:23:00.019538",
                "units": "counts (raw, no instrument correction)",
                "source": "Sigs.mat variable sig1 from the original project; station/time from Kalman_Earth.m",
            },
            do_compression=True,
        )
    print("data written to", ROOT / "data")


if __name__ == "__main__":
    main()
