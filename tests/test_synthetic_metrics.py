import numpy as np

from kalman_denoise.io import load_config
from kalman_denoise.metrics import snr_db, snr_db_original
from kalman_denoise.synthetic import chirp_signal, hyperbolic_gather, ricker


def test_ricker_shape():
    f0 = 20.0
    assert ricker(0.0, f0) == 1.0
    tz = np.sqrt(0.5) / (np.pi * f0)
    assert abs(ricker(tz, f0)) < 1e-15


def test_gather_geometry():
    cfg = load_config()
    d, t, h = hyperbolic_gather(cfg["synthetic"])
    assert d.shape == (1501, 100) and t[-1] == 3.0
    assert np.allclose(d, d[:, ::-1])                    # split spread is symmetric
    j = np.argmin(np.abs(h - 20.0))                      # apex of the first event
    i = np.argmax(d[:350, j])                       # first 0.7 s only
    assert abs(t[i] - np.sqrt(0.5 ** 2 + (20 / 800) ** 2)) <= 0.002


def test_chirp_frequency_law():
    x, t, f = chirp_signal(0.002, 4.0, 1.0, 120.0)
    assert x.size == 2001 and f[0] == 1.0 and abs(f[-1] - 120.0) < 1e-12


def test_snr_definitions():
    c = np.sin(np.linspace(0, 20, 1000))
    assert abs(snr_db(c, 0.5 * c) - 20 * np.log10(2)) < 1e-12
    assert abs(snr_db_original(c, 0.5 * c)) < 1e-12      # the flawed definition gives 0 dB
