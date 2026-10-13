"""Write tests/reference/python_reference.mat, used by the MATLAB parity tests.

    python tools/make_reference.py
"""
import sys
from pathlib import Path

import numpy as np
from scipy.io import savemat

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from kalman_denoise import pipeline as P  # noqa: E402
from kalman_denoise.io import load_config  # noqa: E402
from kalman_denoise.kalman import kalman1d, kalman_rts  # noqa: E402
from kalman_denoise.spectra import amplitude_spectrum, fk_spectrum  # noqa: E402

cfg = load_config()
cols = np.array([0, 49, 99])
ref = {"cols": cols + 1}
s_lo, s_hi = cfg["noise"]["sigma_levels"]
for tag, s in (("lo", s_lo), ("hi", s_hi)):
    r = P.run_synthetic(cfg, s)
    ref[f"snr_{tag}"] = np.array([r["snr"][k] for k in ("input",) + P.METHODS])
    ref[f"snr_orig_{tag}"] = np.array([r["snr_original_definition"][k] for k in ("input",) + P.METHODS])
    if tag == "hi":
        ref["clean_cols"] = r["clean"][:, cols]
        ref["noisy_cols"] = r["noisy"][:, cols]
        for i, m in enumerate(P.METHODS):
            ref[f"out{i + 1}_cols"] = r["out"][m][:, cols]
        f, A = amplitude_spectrum(r["clean"], cfg["synthetic"]["dt_s"])
        ref["amp_f"], ref["amp_A"] = f, A
        k, ff, S = fk_spectrum(r["clean"], cfg["synthetic"]["dt_s"], cfg["synthetic"]["offset_step_m"])
        ref["fk_k"], ref["fk_f"], ref["fk_S_col_k0"] = k, ff, S[:, k.size // 2]
sw = P.run_sweep(cfg)
ref["sweep_Q"] = sw["Q"]
for tag, s in (("lo", s_lo), ("hi", s_hi)):
    ref[f"sweep_fwd_{tag}"] = sw[f"forward_{s:g}"]
    ref[f"sweep_rts_{tag}"] = sw[f"rts_{s:g}"]
ch = P.run_chirp(cfg)
k = cfg["kalman"]
ref["chirp_noisy"] = ch["noisy"]
x, K, Pp = kalman1d(ch["noisy"], k["Q"], k["R"], k["P0"], k["x0"], return_all=True)
ref["chirp_kf"], ref["chirp_K"], ref["chirp_P"] = x, K, Pp
ref["chirp_rts"] = kalman_rts(ch["noisy"], k["Q"], k["R"], k["P0"], k["x0"])
ref["chirp_gain_error"] = ch["max_gain_error"]
eq = P.run_earthquake(cfg)
ref["eq_K"] = np.array([eq["main"]["K"], eq["alt"]["K"]])
ref["eq_f3db"] = np.array([eq["main"]["f3db_hz"], eq["alt"]["f3db_hz"]])
ref["eq_rms_ratio"] = np.array([eq["main"]["rms_ratio"], eq["alt"]["rms_ratio"]])
ref["eq_y_main_head"] = eq["main"]["y"][:2000]
savemat(ROOT / "tests" / "reference" / "python_reference.mat", ref, do_compression=True)
print("wrote", ROOT / "tests" / "reference" / "python_reference.mat")
