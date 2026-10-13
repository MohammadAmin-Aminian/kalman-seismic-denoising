import numpy as np
from scipy.signal import savgol_filter as scipy_savgol

from kalman_denoise.baselines import fx_decon, savgol_filter, svd_denoise
from kalman_denoise.metrics import snr_db
from kalman_denoise.synthetic import ricker


def test_savgol_reproduces_polynomials_exactly():
    t = np.linspace(-1, 1, 200)
    p = 1 - 2 * t + 0.5 * t ** 3 - t ** 4
    assert np.allclose(savgol_filter(p, 4, 27), p, atol=1e-10)


def test_savgol_matches_scipy_interp_mode():
    rng = np.random.default_rng(3)
    x = rng.standard_normal((400, 3))
    ref = scipy_savgol(x, 27, 4, axis=0, mode="interp")
    assert np.allclose(savgol_filter(x, 4, 27), ref, atol=1e-10)


def test_fx_decon_improves_plane_waves():
    dt = 0.002
    t = np.arange(1001) * dt
    x = np.arange(60) * 10.0
    d = sum(a * ricker(t[:, None] - (t0 + p * x[None, :]), 25)
            for a, t0, p in [(1, .3, 4e-4), (-.8, .9, -2e-4), (.6, 1.4, 0.0)])
    noisy = d + 0.2 * np.random.default_rng(0).standard_normal(d.shape)
    out = fx_decon(noisy, dt, 8, 1.0, 0.0, 100.0)
    assert out.shape == d.shape and np.isrealobj(out)
    assert snr_db(d, out) - snr_db(d, noisy) > 8.0


def test_svd_recovers_rank_one():
    u, v = np.arange(1.0, 51.0), np.cos(np.arange(20.0))
    A = np.outer(u, v)
    assert np.allclose(svd_denoise(A, 1), A)
