"""Build the files in data/ (run once; the outputs are committed).

    python tools/make_data.py --sigs Sigs.mat --shot20 shot20.sgy --shot1 shot1.sgy

(each original file is optional; the synthetic and chirp noise are always rebuilt
identically from the fixed seed)

* data/synthetic_noise.mat  unit-variance Gaussian noise (nt x ntraces), single precision
* data/chirp_noise.mat      unit-variance Gaussian noise for the chirp test, single precision
* data/earthquake_IU_DAV_BHZ.mat  vertical broadband record (counts) taken from Sigs.mat/sig1
* data/modelled_shot20.mat  shot 20 of the original modelled file shot20.sgy (noise-free)
* data/modelled_noise.mat   unit-variance Gaussian noise for the modelled shot (single precision)
* data/field_shot_ffid3043.mat  field record FFID 3043 of the original shot1.sgy

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
sys.path.insert(0, str(ROOT / "tools"))
from segy import read_segy  # noqa: E402

SEED = 20210926  # date of the original synthetic experiment script


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sigs", help="original Sigs.mat (only needed for the earthquake file)")
    ap.add_argument("--shot20", help="original shot20.sgy (modelled shots)")
    ap.add_argument("--shot1", help="original shot1.sgy (field record)")
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
    if a.shot20:
        # little-endian file whose binary header says IBM (1) but whose samples are IEEE floats
        d, h, b, _ = read_segy(a.shot20, endian="<", sample_format=5)
        m = h["fldr"] == 20
        sc = 10.0  # coordinate scalar -10 -> divide by 10
        offs = (h["gx"][m] - h["sx"][m]) / sc
        savemat(ROOT / cfg["modelled"]["file"], {
            "d": d[:, m].astype(np.float32), "offset_m": offs[None, :], "dt_s": b["dt_us"] * 1e-6,
            "source_x_m": h["sx"][m][0] / sc,
            "source": "shot20.sgy (original project, folder Data_Syn), shot (FFID) 20 of 20; noise-free modelled data",
        }, do_compression=True)
        mn = np.random.default_rng(SEED + 1).standard_normal((d.shape[0], int(m.sum())))
        savemat(ROOT / cfg["modelled"]["noise_file"], {"noise": mn.astype(np.float32), "seed": SEED + 1},
                do_compression=True)
    if a.shot1:
        d, h, b, text = read_segy(a.shot1, endian=">")
        m = h["fldr"] == 3043
        savemat(ROOT / cfg["field"]["file"], {
            "d": d[:, m].astype(np.float32), "offset_m": h["offset"][m][None, :].astype(float),
            "dt_s": b["dt_us"] * 1e-6, "ffid": 3043,
            "source": "shot1.sgy (original project), field record FFID 3043; area Tang-e-Homam, line ns02",
        }, do_compression=True)
    print("data written to", ROOT / "data")


if __name__ == "__main__":
    main()
