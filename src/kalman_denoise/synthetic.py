"""Deterministic synthetic signals (no random numbers are drawn here)."""
from __future__ import annotations

import numpy as np

__all__ = ["ricker", "time_axis", "offsets", "hyperbolic_gather", "chirp_signal"]


def ricker(t, f0):
    """Ricker wavelet with peak frequency f0 (Hz), unit amplitude at t = 0."""
    a = (np.pi * f0 * np.asarray(t, dtype=float)) ** 2
    return (1.0 - 2.0 * a) * np.exp(-a)


def time_axis(dt, tmax):
    nt = int(np.floor(tmax / dt + 1e-9)) + 1
    return np.arange(nt) * dt


def offsets(cfg_syn):
    """Offsets (m) of the gather; split spread [-h_max..-h_min, h_min..h_max] if mirrored."""
    s = cfg_syn
    n = int(round((s["offset_stop_m"] - s["offset_start_m"]) / s["offset_step_m"])) + 1
    h = s["offset_start_m"] + s["offset_step_m"] * np.arange(n)
    if s.get("mirror", False):
        h = np.concatenate([-h[::-1], h])
    return h


def hyperbolic_gather(cfg_syn):
    """Clean CMP-like gather of hyperbolic reflections.

    d(t, h) = sum_j amp_j * ricker(t - sqrt(tau_j^2 + (h / v_j)^2)).
    Returns (data[nt, nh], t[nt], h[nh]).
    """
    s = cfg_syn
    t = time_axis(s["dt_s"], s["tmax_s"])
    h = offsets(s)
    d = np.zeros((t.size, h.size))
    for tau, v, a in zip(s["tau_s"], s["v_mps"], s["amp"]):
        tt = np.sqrt(tau ** 2 + (h / v) ** 2)
        d += a * ricker(t[:, None] - tt[None, :], s["f0_hz"])
    return d, t, h


def chirp_signal(dt, duration, f_start, f_end):
    """Unit-amplitude linear chirp sin(2*pi*(f_start*t + 0.5*rate*t^2)).

    Returns (signal, t, instantaneous_frequency).
    """
    t = time_axis(dt, duration)
    rate = (f_end - f_start) / duration
    phase = 2.0 * np.pi * (f_start * t + 0.5 * rate * t ** 2)
    return np.sin(phase), t, f_start + rate * t
