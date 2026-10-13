"""Scalar random-walk Kalman filter applied trace by trace.

State-space model (one model per trace, time index k):

    x_k = x_{k-1} + w_k,   w_k ~ N(0, Q)     (state: the clean sample)
    z_k = x_k + v_k,       v_k ~ N(0, R)     (observation: the recorded sample)

Filter recursion (Welch & Bishop, 2006, eqs. 1.9-1.13 with A = H = 1):

    prediction:  x_k^- = x_{k-1}^+,          P_k^- = P_{k-1}^+ + Q
    gain:        K_k   = P_k^- / (P_k^- + R)
    update:      x_k^+ = x_k^- + K_k (z_k - x_k^-),   P_k^+ = (1 - K_k) P_k^-

Missing samples (NaN) are handled as missing observations: only the
prediction step is applied.
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "kalman1d",
    "kalman_rts",
    "steady_state",
    "ewma_response",
    "cutoff_frequency",
]


def _as_2d(z):
    z = np.asarray(z, dtype=float)
    was_1d = z.ndim == 1
    if was_1d:
        z = z[:, None]
    if z.ndim != 2:
        raise ValueError("input must be a vector or a 2-D array (time x traces)")
    return z, was_1d


def _check(Q, R, P0):
    if not (Q >= 0 and R > 0 and P0 >= 0):
        raise ValueError("require Q >= 0, R > 0 and P0 >= 0")


def kalman1d(z, Q, R=0.05, P0=1.01, x0=0.0, return_all=False):
    """Forward (causal) random-walk Kalman filter along axis 0.

    Parameters
    ----------
    z : array_like, shape (nt,) or (nt, ntraces)
        Noisy observations; each column is filtered independently.
    Q, R : float
        Process- and measurement-noise variances.
    P0, x0 : float
        Prior variance and prior mean of the state before the first sample.
    return_all : bool
        Also return the gain K and posterior variance P (same shape as z).

    Returns
    -------
    x : ndarray, same shape as z
        Filtered estimate x_k^+.
    """
    _check(Q, R, P0)
    z2, was_1d = _as_2d(z)
    nt, ntr = z2.shape
    x = np.empty_like(z2)
    K = np.zeros_like(z2)
    P = np.empty_like(z2)
    xk = np.full(ntr, float(x0))
    Pk = np.full(ntr, float(P0))
    for k in range(nt):
        Pm = Pk + Q if k > 0 else Pk.copy()  # the prior already describes sample 0
        obs = ~np.isnan(z2[k])
        Kk = np.where(obs, Pm / (Pm + R), 0.0)
        innov = np.where(obs, z2[k] - xk, 0.0)
        xk = xk + Kk * innov
        Pk = (1.0 - Kk) * Pm
        x[k], K[k], P[k] = xk, Kk, Pk
    if was_1d:
        x, K, P = x[:, 0], K[:, 0], P[:, 0]
    return (x, K, P) if return_all else x


def kalman_rts(z, Q, R=0.05, P0=1.01, x0=0.0):
    """Rauch-Tung-Striebel (fixed-interval) smoother for the same model.

    Extension, not part of the original study: runs the forward filter and a
    backward pass, which removes the phase delay of the causal filter.
    """
    xf, _, Pf = kalman1d(z, Q, R, P0, x0, return_all=True)
    xf2, was_1d = _as_2d(xf)
    Pf2, _ = _as_2d(Pf)
    xs = xf2.copy()
    for k in range(xs.shape[0] - 2, -1, -1):
        Pp = Pf2[k] + Q                       # P_{k+1}^-
        C = np.divide(Pf2[k], Pp, out=np.zeros_like(Pp), where=Pp > 0)
        xs[k] = xf2[k] + C * (xs[k + 1] - xf2[k])
    return xs[:, 0] if was_1d else xs


def steady_state(Q, R=0.05):
    """Steady-state prior variance, posterior variance and gain.

    P^- solves P^-^2 - Q P^- - Q R = 0, so
    P^- = (Q + sqrt(Q^2 + 4 Q R)) / 2,  K = P^- / (P^- + R),  P^+ = (1 - K) P^-.
    """
    _check(Q, R, 0.0)
    Pm = 0.5 * (Q + np.sqrt(Q * Q + 4.0 * Q * R))
    K = Pm / (Pm + R)
    return {"P_prior": Pm, "P_post": (1.0 - K) * Pm, "K": K}


def ewma_response(K, f, dt):
    """Complex frequency response of the steady-state filter.

    In steady state the update is x_k = (1 - K) x_{k-1} + K z_k, i.e.
    H(z) = K / (1 - (1 - K) z^-1).
    """
    w = 2.0 * np.pi * np.asarray(f, dtype=float) * dt
    return K / (1.0 - (1.0 - K) * np.exp(-1j * w))


def cutoff_frequency(K, dt):
    """-3 dB frequency (Hz) of the steady-state filter (nan if none below Nyquist)."""
    a = 1.0 - K
    if a <= 0:
        return np.nan
    c = (1.0 + a * a - 2.0 * K * K) / (2.0 * a)
    if c < -1.0:
        return np.nan
    return np.arccos(min(c, 1.0)) / (2.0 * np.pi * dt)
