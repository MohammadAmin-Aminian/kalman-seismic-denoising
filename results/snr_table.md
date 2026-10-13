Output SNR (dB) = 10 log10(Σ clean² / Σ (output − clean)²), synthetic gather.

| Method | σ = 0.04 | σ = 0.4 |
|---|---:|---:|
| Input (noisy) | 21.48 | 1.48 |
| Savitzky-Golay | 15.69 | 9.16 |
| f-x decon | 4.85 | 3.89 |
| SVD | 7.67 | 5.05 |
| Kalman | 4.67 | 3.56 |
| Kalman RTS | 6.45 | 5.58 |
| Kalman, best Q (sweep) | 21.78 | 5.23 |
| Kalman RTS, best Q (sweep) | 23.75 | 8.32 |
