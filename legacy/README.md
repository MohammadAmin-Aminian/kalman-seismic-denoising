# Original scripts (2021) — for reference only

These are the original MATLAB scripts of the study, kept so that every
correction in [`docs/CORRECTIONS.md`](../docs/CORRECTIONS.md) can be traced.
They are **not** maintained and **do not run** as they are: they need
third-party toolboxes, files on the author's former computer and functions
that are missing (`kalman`, `VRE`). The only change made to them is that
absolute local file paths were replaced with `<local path>`.

| File | Purpose |
|---|---|
| `kalman_26_sep.m` | Final synthetic comparison (S-G, f-x, SVD, Kalman). A copy with σ = 0.4 instead of 0.04 existed as `kalman_3.m`. |
| `random_noise_attenuation.m`, `random_noise_attenuation_31_Aug.m` | Earlier versions of the synthetic comparison |
| `kalman_22_sep.m` | Chirp tests |
| `Kalman_Real.m`, `kalman_Real_22_sep.m` | Real shot gathers (SEG-Y files not available) |
| `Kalman_Earth.m`, `kalmancode.m` | Earthquake record |

Functions the scripts call that came from elsewhere, and are therefore not
redistributed here: `Kalman1D.m` (Y. Khmou, MATLAB File Exchange, 2012, after
Welch & Bishop), `hyperbolic_events.m`, `fx_decon.m` and `fk_spectra.m`
(SeismicLab, M. D. Sacchi), `wiggle.m` (R. S. Portugal), `agc.m` and
`altreadsegy.m` (CREWES). The repository re-implements what was needed.
