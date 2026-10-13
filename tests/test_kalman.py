import numpy as np

from kalman_denoise.kalman import (cutoff_frequency, ewma_response, kalman1d, kalman_rts,
                                   steady_state)


def test_gain_converges_to_closed_form():
    for Q in (1e-5, 4e-3, 0.5):
        _, K, _ = kalman1d(np.zeros(4000), Q, 0.05, return_all=True)
        assert abs(K[-1] - steady_state(Q, 0.05)["K"]) < 1e-12


def test_constant_signal_is_recovered():
    x = kalman1d(np.full(500, 3.7), 1e-3, 0.05)
    assert abs(x[-1] - 3.7) < 1e-6


def test_steady_state_is_first_order_recursion():
    rng = np.random.default_rng(0)
    z = rng.standard_normal(3000)
    Q, R = 4e-3, 0.05
    x = kalman1d(z, Q, R)
    K = steady_state(Q, R)["K"]
    y = x[999]
    for k in range(1000, 3000):
        y = (1 - K) * y + K * z[k]
    assert abs(y - x[-1]) < 1e-9


def test_mse_matches_posterior_variance_on_true_model():
    rng = np.random.default_rng(1)
    Q, R, n, m = 0.01, 0.05, 4000, 200
    truth = np.cumsum(np.sqrt(Q) * rng.standard_normal((n, m)), axis=0)
    z = truth + np.sqrt(R) * rng.standard_normal((n, m))
    x = kalman1d(z, Q, R, P0=1e6)
    mse = np.mean((x[500:] - truth[500:]) ** 2)
    P = steady_state(Q, R)["P_post"]
    assert abs(mse / P - 1) < 0.03
    xs = kalman_rts(z, Q, R, P0=1e6)
    assert np.mean((xs[500:-500] - truth[500:-500]) ** 2) < 0.8 * mse


def test_traces_are_independent():
    rng = np.random.default_rng(2)
    z = rng.standard_normal((300, 5))
    X = kalman1d(z, 4e-3)
    for j in range(5):
        assert np.allclose(X[:, j], kalman1d(z[:, j], 4e-3))


def test_missing_samples_are_skipped():
    z = np.sin(np.linspace(0, 10, 400))
    z[100:110] = np.nan
    x, K, _ = kalman1d(z, 4e-3, return_all=True)
    assert np.all(np.isfinite(x)) and np.all(K[100:110] == 0)
    assert np.all(x[100:110] == x[99])


def test_cutoff_frequency_is_minus_3db():
    K, dt = steady_state(4e-3, 0.05)["K"], 0.002
    fc = cutoff_frequency(K, dt)
    assert abs(abs(ewma_response(K, fc, dt)) - 1 / np.sqrt(2)) < 1e-12
    assert abs(fc - 22.5836) < 1e-3


def test_invalid_parameters_raise():
    for bad in ((-1, 0.05), (1e-3, 0.0)):
        try:
            kalman1d(np.zeros(10), *bad)
        except ValueError:
            continue
        raise AssertionError("expected ValueError")
