"""Pre-processing of shot gathers."""
from __future__ import annotations

import numpy as np

__all__ = ["agc", "agc_window_samples"]


def agc_window_samples(window_s, dt):
    """Odd window length (samples) closest to window_s."""
    return 2 * int(round(window_s / dt / 2.0)) + 1


def agc(d, window):
    """Automatic gain control: divide each sample by the RMS in a centred window.

    The window (odd number of samples) shrinks at the trace ends. The gather
    is first scaled by its maximum absolute value (AGC is scale invariant) so
    that the running sums are well conditioned; a stabilisation of 1e-6 times
    the largest RMS avoids division by zero in dead zones.
    """
    d = np.asarray(d, dtype=float)
    one_d = d.ndim == 1
    X = d[:, None] if one_d else d
    if window % 2 == 0 or window < 1:
        raise ValueError("window must be a positive odd number of samples")
    peak = np.max(np.abs(X))
    if peak == 0:
        return d.copy()
    X = X / peak
    nt = X.shape[0]
    h = window // 2
    c = np.vstack([np.zeros((1, X.shape[1])), np.cumsum(X ** 2, axis=0)])
    i = np.arange(nt)
    lo = np.maximum(i - h, 0)
    hi = np.minimum(i + h, nt - 1)
    ms = (c[hi + 1] - c[lo]) / (hi - lo + 1)[:, None]
    rms = np.sqrt(np.maximum(ms, 0.0))
    y = X / (rms + 1e-6 * np.max(rms))
    return y[:, 0] if one_d else y
