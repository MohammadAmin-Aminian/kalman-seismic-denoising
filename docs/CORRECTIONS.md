# Corrections relative to the original code

There is no paper for this project, so every entry below is relative to the
original MATLAB scripts in [`legacy/`](../legacy) (2021). Each item lists
what was wrong, what changed and how much it moves the results. The
size of each effect is measured with
[`tools/original_port.py`](../tools/original_port.py), a faithful Python port
of the original code that keeps its bugs
([`results/original_comparison.md`](../results/original_comparison.md)).

## A. Errors that change reported numbers

| # | Original code | Problem | Correction | Effect |
|---|---|---|---|---|
| A1 | `kalman_26_sep.m`, `random_noise_attenuation*.m`, `kalman_22_sep.m` | The SNR was `10 log10(Σ out² / Σ (out − clean)²)`, with the **output** in the numerator instead of the clean signal. A filter that adds energy scores higher, and an output at half amplitude scores 0 dB. | `snr_db = 10 log10(Σ clean² / Σ (out − clean)²)`. The old formula is kept as `snr_db_original` for comparison only. | Kalman at σ = 0.04: 2.08 → 4.67 dB. f-x at σ = 0.4: 0.26 → 3.89 dB. All methods change. |
| A2 | `kalman_Real_22_sep.m` | The "SNR" of the real data was computed against the **noisy input** (`res = xf − shot`). Without a clean reference, this measures how much the filter changed the data, not how good the result is. | Removed. On field data the repository reports filter characteristics and spectra instead. | Number no longer reported. |
| A3 | `kalman_26_sep.m`, `random_noise_attenuation_31_Aug.m`, `kalman_*_sep.m` | The frequency axis of the spectra used `dt1 = 0.001` s and a bin spacing of 1 Hz, but the data are sampled at 2 ms with 1501 samples (true spacing 0.333 Hz). | `amplitude_spectrum` uses `f = k / (N dt)`. | Original frequency axes are scaled by 3: the 20 Hz Ricker peak appeared at 60 Hz. |
| A4 | `kalman_26_sep.m`, `Kalman_Real.m` | The f-k spectrum used `dt = 0.001` s and `dx = 12.5` m. The synthetic gather has 2 ms and 20 m. | `fk_spectrum(d, dt, dx)` takes both from the config (2 ms, 20 m). | Original f axis doubled, k axis scaled by 1.6. |
| A5 | `kalman_26_sep.m` | `fx_decon(data, .001, 8, 0, 0, 200)` was called with the wrong `dt` (1 ms for 2 ms data), so the requested 0–200 Hz band was actually 0–100 Hz. | `dt` is now correct and the band is set explicitly to `fmax_hz = 100`, which keeps the band the original actually used. | None (identical output), but the parameter now means what it says. |
| A6 | `kalmancode.m` | Calls `kalman(sig1, 5e-5)`, which doesn't exist in the project. | Interpreted as `Kalman1D` with Q = 5·10⁻⁵ (confirmed default). | Earthquake figure reproducible. |
| A7 | `Kalman_Earth.m`, `kalmancode.m` | Spectra of the 40 Hz earthquake record were plotted with `pspectrum` normalised frequency × 1000, i.e. 0–1000 "Hz" for data whose Nyquist frequency is 20 Hz. `kalmancode.m` plotted the complex FFT, which shows only its real part. | One-sided amplitude spectrum, 0–20 Hz, log-log axes. | Earthquake spectra now physically meaningful. |
| A8 | `Kalman_Earth.m`, `kalmancode.m` | f-k spectra of a **single trace** were computed (`fk_spectra(Data', …)`). A wavenumber needs several traces, so the result is meaningless. | Removed. | — |
| A9 | `Kalman_Earth.m` | Axis labels said "Amplitude (m)" and "Time (s)", but the data are raw counts plotted against sample number. | Counts against time in seconds (40 Hz). | Labels only. |

## B. Bugs and robustness fixes in the method

| # | Original | Problem | Correction | Effect |
|---|---|---|---|---|
| B1 | `Kalman1D.m` | `N = length(X)` was computed **before** NaN/Inf samples were removed. Any NaN shortens the trace, so the loop reads past its end and every later sample shifts in time. | NaN samples are treated as missing observations: prediction only, gain 0. | None on clean inputs. Correct behaviour with gaps. |
| B2 | `Kalman1D.m` | The first sample was never updated: `Y(1) = 0` always. | The measurement update is also applied to sample 1, using the prior (x0, P0). | ≤ 0.01 dB. |
| B3 | `Kalman_Earth.m` | The record has a DC offset of about 3,900 counts and the filter starts at x0 = 0, which causes a start-up transient. | The record is demeaned first. | Removes the transient. |
| B4 | `hyperbolic_events.m` (SeismicLab) | The "clean" gather contained extra band-limited random noise (`snr = 1000, L = 1000`, unseeded). It was also built in the frequency domain with a truncated wavelet. | The clean gather is now an exact sum of analytic Ricker wavelets: deterministic, and identical in both languages. | < 0.01 dB (Part A of the comparison). |
| B5 | all scripts | `randn` was never seeded, so no run could be repeated exactly. | One noise field is drawn with a fixed seed and **stored** in `data/*.mat`. Python and MATLAB read the same numbers. | Exact reproducibility and cross-language parity. |
| B6 | `kalman_22_sep.m` (chirp) | `chirp(t, f0, t1, f1)` used `f0` left over in the workspace, and `t1 = 0.05` s. The sweep therefore continued linearly to about 320 Hz while sampling at 100 Hz, so it was aliased. The noise was `0.8*rand`, which is uniform with mean 0.4, but labelled "Gaussian". | A properly sampled 1–120 Hz linear chirp at 500 Hz plus zero-mean Gaussian noise. It now serves as a known-answer test of the filter's frequency response. | New validation (Fig. 9). |
| B7 | `kalman_Real_22_sep.m` | Uses `D1`, `D` and `dtsec` before they are defined, so the script can't run. | Rewritten as `run_field`: AGC, then all methods on field record FFID 3043 of `shot1.sgy`. | Fig. 12–13. |
| B10 | `Kalman_Real.m`, `kalman_Real_22_sep.m` | The f-x and f-k calls used `dt = 0.001` s for the 2 ms modelled shots and the 4 ms field record. They also used `dx = 12.5` m, while the modelled shots are 20 m and the field record 25 m. | Each dataset uses its own dt (from the SEG-Y header) and trace spacing. | Correct frequency bands. |
| B11 | `shot20.sgy` | The modelled file is little-endian and its binary header declares format 1 (IBM float), but the samples are IEEE floats. Read as IBM, the values are meaningless (±3.6·10¹⁶). | `tools/segy.py` reads it as little-endian IEEE. | Data decoded correctly. |
| B12 | `Kalman_Real.m` | The modelled shots are noise-free, so "denoising" them measured nothing. | Gaussian noise (0 dB, stored) is added to shot 20 after AGC, and every method is scored against the noise-free gather. | New benchmark (Fig. 11). |
| B8 | `kalman_26_sep.m`, `random_noise_attenuation*.m` | SVD denoising kept 18 or 50 full-size rank-one images in memory (`Image(:,:,i)`) and summed them. | Direct truncated SVD `U(:,1:r) S V(:,1:r)'`. | Identical result, far less memory. |
| B9 | `kalman_26_sep.m` | `VRE(xf, 2, 10)` ("Kalman + VRE") calls a function that isn't in the project. | Dropped (confirmed). | — |

## C. Toolbox and portability changes

| # | Original | Replacement (base MATLAB / numpy-scipy) | Check |
|---|---|---|---|
| C1 | `sgolayfilt` (Signal Processing Toolbox) | `savgol_smooth` / `savgol_filter`, with the same edge rule | Equal to `scipy.signal.savgol_filter(mode="interp")` to 1e-10. Reproduces degree-4 polynomials exactly. |
| C2 | `fx_decon` (SeismicLab) | `fx_decon_filter` / `fx_decon`, a new implementation of Canales (1984) | Same output as the original port. +11 dB on plane waves (unit test). |
| C3 | `fk_spectra`, `wiggle`, `agc`, `altreadsegy`, `ricker` (third-party) | `fk_spectrum`, image displays, analytic `ricker_wavelet` | Parity tests between Python and MATLAB. |
| C4 | `chirp`, `hilbert`, `pspectrum`, `normalize` (Signal Processing Toolbox) | `chirp_signal`, `analytic_signal` (FFT), `amplitude_spectrum` | Unit tests. |
| C5 | Hard-coded Windows paths (`C:\Users\…`) | Paths relative to the repository root and one JSON config | Runs anywhere. |
| C6 | `agc(data,100,1,0)`, `altreadsegy` (third-party) | `agc` / `agc_gain`: centred RMS AGC, 100-sample window, shrinking at trace ends. `tools/segy.py` handles the one-off SEG-Y conversion. | Unit tests (RMS normalisation, scale invariance). The AGC definition is an assumption. |

## D. Choices that were ambiguous in the originals

| # | Ambiguity | Decision |
|---|---|---|
| D1 | Noise level: `kalman_3.m` uses σ = 0.4, the later `kalman_26_sep.m` uses σ = 0.04 | Both are reported. σ = 0.04 is the "main" level in the config. |
| D2 | Earthquake Q: `kalmancode.m` uses 5·10⁻⁵, `Kalman_Earth.m` uses 0.005 | Both are shown (`earthquake.Q`, `earthquake.Q_alt`). |
| D3 | Spectra used the **maximum** over traces | The **mean** over traces, which is less sensitive to single noisy traces. Display only. |
| D4 | `Sigs.mat` holds `sig1` (int32) and `sig2` (int16). They correlate at r = 0.986 and `sig2` is unidentified. | Only `sig1` is published. Its station and start time come from the file name in `Kalman_Earth.m`. |
| D5 | SNR values written in the titles of `random_noise_attenuation.m` (noisy 4, S-G 13, FX 10, SVD 6, Kalman 7 dB) | Only Kalman (6.86 dB) can be reproduced from that script. The others don't match its parameters and were probably typed in from other runs. Reported as found. |

## E. Additions (not in the original study, clearly labelled)

* **Kalman RTS smoother** (`kalman_rts`): the fixed-interval smoother for the same model. It removes the phase lag of the causal filter. It appears as a separate column in the results and never replaces the original method.
* **Q sweep** (Fig. 7), **steady-state analysis** (Fig. 8) and the **chirp known-answer test** (Fig. 9).
