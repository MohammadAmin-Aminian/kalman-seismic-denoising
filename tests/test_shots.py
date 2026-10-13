import json

import numpy as np

from kalman_denoise import pipeline as P
from kalman_denoise.io import load_config, repo_root
from kalman_denoise.preprocessing import agc, agc_window_samples


def test_agc_normalises_rms():
    t = np.arange(4000) * 0.002
    x = np.exp(-t) * np.sin(2 * np.pi * 15 * t) * 1e6          # decaying, large amplitude
    y = agc(x, agc_window_samples(0.2, 0.002))
    rms = np.sqrt(np.convolve(y ** 2, np.ones(101) / 101, mode="valid"))
    assert np.allclose(rms, 1.0, atol=0.05)


def test_agc_is_scale_invariant_and_handles_zeros():
    rng = np.random.default_rng(0)
    d = rng.standard_normal((500, 4))
    d[:100] = 0.0
    a, b = agc(d, 51), agc(1e9 * d, 51)
    assert np.allclose(a, b) and np.all(a[:70] == 0)


def test_agc_window_must_be_odd():
    try:
        agc(np.ones(10), 4)
    except ValueError:
        return
    raise AssertionError("expected ValueError")


def test_shot_numbers_are_reproduced():
    cfg = load_config()
    with open(repo_root() / "results" / "metrics.json", encoding="utf-8") as fh:
        pub = json.load(fh)
    mo = P.run_modelled(cfg)
    for k, v in pub["modelled_shot"]["snr_db"].items():
        assert abs(mo["snr"][k] - v) < 1e-6, k
    fd = P.run_field(cfg)
    for k, v in pub["field_record"]["removed_energy_fraction"].items():
        assert abs(fd["removed_fraction"][k] - v) < 1e-9, k
