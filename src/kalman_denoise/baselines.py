"""Reference denoising methods used for comparison with the Kalman filter."""
from __future__ import annotations

import numpy as np

__all__ = ["savgol_coefficients", "savgol_filter", "fx_decon", "svd_denoise"]


def savgol_coefficients(order, frame):
    """Projection matrix B (frame x frame) of a least-squares polynomial fit.

    Row i gives the weights that evaluate, at window position i, the
    polynomial of degree `order` fitted to the whole window.
    """
    if frame % 2 == 0 or frame <= order:
        raise ValueError("frame must be odd and larger than order")
    m = (frame - 1) // 2
    s = np.arange(-m, m + 1) / m                    # scaled positions for conditioning
    V = np.vander(s, order + 1, increasing=True)
    return V @ np.linalg.solve(V.T @ V, V.T)


def savgol_filter(x, order, frame):
    """Savitzky-Golay smoothing along axis 0 (same edge rule as MATLAB sgolayfilt).

    Interior samples use the central row of B; the first/last (frame-1)/2
    samples are evaluated from the polynomial fitted to the first/last frame.
    """
    x = np.asarray(x, dtype=float)
    one_d = x.ndim == 1
    X = x[:, None] if one_d else x
    nt = X.shape[0]
    if nt < frame:
        raise ValueError("signal shorter than the frame")
    B = savgol_coefficients(order, frame)
    m = (frame - 1) // 2
    c = B[m]
    y = np.empty_like(X)
    # interior: correlation with the central row
    from numpy.lib.stride_tricks import sliding_window_view
    win = sliding_window_view(X, frame, axis=0)     # (nt-frame+1, ntr, frame)
    y[m:nt - m] = win @ c
    y[:m] = B[:m] @ X[:frame]
    y[nt - m:] = B[m + 1:] @ X[nt - frame:]
    return y[:, 0] if one_d else y


def _ar_predict(x, lf, mu):
    """Forward and backward one-step linear prediction of a complex spatial series."""
    n = x.size
    yf = np.zeros(n, dtype=complex)
    yb = np.zeros(n, dtype=complex)
    if n <= lf:
        return yf, yb
    # forward: x[i] ~ sum_j a_j x[i-j], j = 1..lf, for i = lf..n-1
    Mf = np.column_stack([x[lf - j:n - j] for j in range(1, lf + 1)])
    # backward: x[i] ~ sum_j b_j x[i+j], j = 1..lf, for i = 0..n-lf-1
    Mb = np.column_stack([x[j:n - lf + j] for j in range(1, lf + 1)])
    for M, target, sl, out in ((Mf, x[lf:], slice(lf, n), yf), (Mb, x[:n - lf], slice(0, n - lf), yb)):
        A = M.conj().T @ M
        beta = np.real(A[0, 0]) * mu / 100.0
        a = np.linalg.solve(A + beta * np.eye(lf), M.conj().T @ target)
        out[sl] = M @ a
    return yf, yb


def fx_decon(d, dt, lf=8, mu=0.0, fmin=0.0, fmax=None):
    """f-x deconvolution (Canales, 1984): spatial prediction filtering per frequency.

    For every frequency in [fmin, fmax] a forward and a backward prediction
    filter of length `lf` are estimated by damped least squares (damping
    mu % of the zero-lag autocorrelation) and the predictions are averaged.
    Frequencies outside [fmin, fmax] are set to zero, as in SeismicLab.
    """
    d = np.asarray(d, dtype=float)
    nt, ntr = d.shape
    nf = 1 << int(np.ceil(np.log2(nt)))
    if fmax is None:
        fmax = 0.5 / dt
    ilow = max(int(np.floor(fmin * dt * nf)), 0)
    ihigh = min(int(np.floor(fmax * dt * nf)), nf // 2)
    D = np.fft.fft(d, nf, axis=0)
    F = np.zeros_like(D)
    Bk = np.zeros_like(D)
    for k in range(ilow, ihigh + 1):
        F[k], Bk[k] = _ar_predict(D[k], lf, mu)
    for k in range(nf // 2 + 1, nf):             # Hermitian symmetry -> real output
        F[k] = np.conj(F[nf - k])
        Bk[k] = np.conj(Bk[nf - k])
    out = np.real(np.fft.ifft(F, axis=0))[:nt] + np.real(np.fft.ifft(Bk, axis=0))[:nt]
    out[:, lf:ntr - lf] /= 2.0                  # interior traces have both predictions
    return out


def svd_denoise(d, rank):
    """Rank-r reconstruction (sum of the first r eigenimages)."""
    U, s, Vt = np.linalg.svd(np.asarray(d, dtype=float), full_matrices=False)
    return (U[:, :rank] * s[:rank]) @ Vt[:rank]
