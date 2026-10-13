"""Spectral displays."""
from __future__ import annotations

import numpy as np

__all__ = ["amplitude_spectrum", "fk_spectrum"]


def amplitude_spectrum(d, dt):
    """One-sided amplitude spectrum |FFT| averaged over traces.

    Returns (f[Hz], A) with f = 0 .. Nyquist.
    """
    d = np.asarray(d, dtype=float)
    if d.ndim == 1:
        d = d[:, None]
    nt = d.shape[0]
    A = np.abs(np.fft.rfft(d, axis=0)).mean(axis=1)
    return np.fft.rfftfreq(nt, dt), A


def fk_spectrum(d, dt, dx):
    """Centred f-k amplitude spectrum (dB, relative to its maximum).

    Returns (k[1/km], f[Hz] >= 0, S[nf, nk]).
    """
    d = np.asarray(d, dtype=float)
    nt, nx = d.shape
    nfft_t = 1 << int(np.ceil(np.log2(nt)))
    nfft_x = 1 << int(np.ceil(np.log2(nx)))
    D = np.fft.fftshift(np.fft.fft2(d, (nfft_t, nfft_x)), axes=1)
    D = np.abs(D[: nfft_t // 2 + 1])
    S = 20.0 * np.log10(D / D.max() + 1e-12)
    f = np.arange(nfft_t // 2 + 1) / (nfft_t * dt)
    k = (np.arange(nfft_x) - nfft_x // 2) / (nfft_x * dx) * 1000.0
    return k, f, S
