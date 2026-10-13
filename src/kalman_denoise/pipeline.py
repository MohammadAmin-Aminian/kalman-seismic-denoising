"""End-to-end experiments. Each function returns plain dicts of arrays/numbers."""
from __future__ import annotations

import numpy as np

from .baselines import fx_decon, savgol_filter, svd_denoise
from .io import load_earthquake, load_noise, load_shot
from .kalman import cutoff_frequency, ewma_response, kalman1d, kalman_rts, steady_state
from .metrics import snr_db, snr_db_original
from .preprocessing import agc, agc_window_samples
from .synthetic import chirp_signal, hyperbolic_gather

METHODS = ("Savitzky-Golay", "f-x decon", "SVD", "Kalman", "Kalman RTS")


def analytic_signal(x):
    """Analytic signal via FFT (same as scipy.signal.hilbert)."""
    x = np.asarray(x, dtype=float)
    n = x.size
    X = np.fft.fft(x)
    h = np.zeros(n)
    h[0] = 1.0
    if n % 2 == 0:
        h[n // 2] = 1.0
        h[1:n // 2] = 2.0
    else:
        h[1:(n + 1) // 2] = 2.0
    return np.fft.ifft(X * h)


def denoise_all(noisy, cfg, Q=None, dt=None, fx=None):
    """Apply every method with the configured parameters (dt, Q and f-x settings can be overridden)."""
    k = cfg["kalman"]
    Q = k["Q"] if Q is None else Q
    dt = cfg["synthetic"]["dt_s"] if dt is None else dt
    fx = cfg["fx"] if fx is None else fx
    return {
        "Savitzky-Golay": savgol_filter(noisy, cfg["savgol"]["order"], cfg["savgol"]["frame"]),
        "f-x decon": fx_decon(noisy, dt, fx["operator_length"], fx["mu_percent"],
                              fx["fmin_hz"], fx["fmax_hz"]),
        "SVD": svd_denoise(noisy, cfg["svd"]["rank"]),
        "Kalman": kalman1d(noisy, Q, k["R"], k["P0"], k["x0"]),
        "Kalman RTS": kalman_rts(noisy, Q, k["R"], k["P0"], k["x0"]),
    }


def run_synthetic(cfg, sigma):
    clean, t, h = hyperbolic_gather(cfg["synthetic"])
    noise = sigma * load_noise(cfg["noise"]["file"])
    noisy = clean + noise
    out = denoise_all(noisy, cfg)
    snr = {"input": snr_db(clean, noisy)}
    snr_orig = {"input": snr_db(clean, noisy)}
    for m, e in out.items():
        snr[m] = snr_db(clean, e)
        snr_orig[m] = snr_db_original(clean, e)
    return {"sigma": sigma, "t": t, "h": h, "clean": clean, "noisy": noisy,
            "out": out, "snr": snr, "snr_original_definition": snr_orig}


def run_sweep(cfg):
    sw, k = cfg["sweep"], cfg["kalman"]
    Q = np.logspace(sw["log10_Q_min"], sw["log10_Q_max"], int(sw["n"]))
    clean, _, _ = hyperbolic_gather(cfg["synthetic"])
    unit = load_noise(cfg["noise"]["file"])
    res = {"Q": Q, "K": np.array([steady_state(q, k["R"])["K"] for q in Q])}
    for sigma in cfg["noise"]["sigma_levels"]:
        noisy = clean + sigma * unit
        fwd = [snr_db(clean, kalman1d(noisy, q, k["R"], k["P0"], k["x0"])) for q in Q]
        rts = [snr_db(clean, kalman_rts(noisy, q, k["R"], k["P0"], k["x0"])) for q in Q]
        res[f"forward_{sigma:g}"] = np.array(fwd)
        res[f"rts_{sigma:g}"] = np.array(rts)
        res[f"input_{sigma:g}"] = snr_db(clean, noisy)
    return res


def run_chirp(cfg):
    c, k = cfg["chirp"], cfg["kalman"]
    x, t, finst = chirp_signal(c["dt_s"], c["duration_s"], c["f_start_hz"], c["f_end_hz"])
    noisy = x + c["sigma"] * load_noise(c["noise_file"])[:, 0]
    y_clean = kalman1d(x, k["Q"], k["R"], k["P0"], k["x0"])
    y_noisy = kalman1d(noisy, k["Q"], k["R"], k["P0"], k["x0"])
    ss = steady_state(k["Q"], k["R"])
    gain_meas = np.abs(analytic_signal(y_clean)) / np.abs(analytic_signal(x))
    gain_theory = np.abs(ewma_response(ss["K"], finst, c["dt_s"]))
    # compare away from both ends (start-up transient and FFT edge effects)
    n = t.size
    mid = slice(int(0.1 * n), int(0.9 * n))
    err = np.max(np.abs(gain_meas[mid] - gain_theory[mid]))
    return {"t": t, "f_inst": finst, "clean": x, "noisy": noisy, "filtered_clean": y_clean,
            "filtered_noisy": y_noisy, "gain_measured": gain_meas, "gain_theory": gain_theory,
            "max_gain_error": float(err), "K": ss["K"],
            "f3db_hz": cutoff_frequency(ss["K"], c["dt_s"])}


def run_earthquake(cfg):
    e = cfg["earthquake"]
    eq = load_earthquake(e["file"])
    x = eq["counts"] - eq["counts"].mean()
    dt = 1.0 / eq["fs_hz"]
    res = {"meta": eq, "dt": dt, "x": x, "t": np.arange(x.size) * dt}
    for name, Q in (("main", e["Q"]), ("alt", e["Q_alt"])):
        y = kalman1d(x, Q, e["R"], cfg["kalman"]["P0"], 0.0)
        ss = steady_state(Q, e["R"])
        res[name] = {"Q": Q, "y": y, "K": ss["K"], "f3db_hz": cutoff_frequency(ss["K"], dt),
                     "rms_ratio": float(np.sqrt(np.mean(y ** 2) / np.mean(x ** 2)))}
    return res


def run_modelled(cfg):
    """Modelled shot gather with ground roll: AGC, add Gaussian noise, denoise, score against the truth."""
    m, k = cfg["modelled"], cfg["kalman"]
    shot = load_shot(m["file"])
    dt = shot["dt_s"]
    clean = agc(shot["d"], agc_window_samples(m["agc_window_s"], dt))
    rms = float(np.sqrt(np.mean(clean ** 2)))
    noisy = clean + m["sigma"] * rms * load_noise(m["noise_file"])
    out = denoise_all(noisy, cfg, Q=m["Q"], dt=dt, fx=m["fx"])
    snr = {"input": snr_db(clean, noisy)}
    snr.update({name: snr_db(clean, e) for name, e in out.items()})
    sw = cfg["sweep"]
    Qs = np.logspace(sw["log10_Q_min"], sw["log10_Q_max"], int(sw["n"]))
    fwd = np.array([snr_db(clean, kalman1d(noisy, q, k["R"], k["P0"], k["x0"])) for q in Qs])
    rts = np.array([snr_db(clean, kalman_rts(noisy, q, k["R"], k["P0"], k["x0"])) for q in Qs])
    return {"t": np.arange(clean.shape[0]) * dt, "h": shot["offset_m"], "dt": dt, "clean": clean,
            "noisy": noisy, "out": out, "snr": snr, "sweep_Q": Qs, "sweep_forward": fwd, "sweep_rts": rts,
            "Q": m["Q"], "sigma": m["sigma"]}


def run_field(cfg):
    """Field record: AGC and denoising; no ground truth, so only removed energy is reported."""
    f, k = cfg["field"], cfg["kalman"]
    shot = load_shot(f["file"])
    dt = shot["dt_s"]
    x = agc(shot["d"], agc_window_samples(f["agc_window_s"], dt))
    out = denoise_all(x, cfg, Q=f["Q"], dt=dt, fx=f["fx"])
    removed = {name: float(np.sum((x - e) ** 2) / np.sum(x ** 2)) for name, e in out.items()}
    ss = steady_state(f["Q"], k["R"])
    return {"t": np.arange(x.shape[0]) * dt, "h": shot["offset_m"], "dt": dt, "x": x, "out": out,
            "removed_fraction": removed, "K": ss["K"], "f3db_hz": cutoff_frequency(ss["K"], dt), "Q": f["Q"]}
