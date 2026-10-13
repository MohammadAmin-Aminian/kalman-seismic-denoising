Output SNR (dB) = 10 log10(Σ clean² / Σ (output − clean)²).

## Synthetic gather

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

## Modelled shot gather (noise σ = 1 × RMS, Kalman Q = 0.005)

| Method | Output SNR (dB) |
|---|---:|
| Input (noisy) | -0.02 |
| Savitzky-Golay | 3.76 |
| f-x decon | 1.20 |
| SVD | 1.35 |
| Kalman | 3.12 |
| Kalman RTS | 4.28 |
| Kalman, best Q (sweep) | 3.60 |
| Kalman RTS, best Q (sweep) | 5.31 |

## Field record FFID 3043 (no ground truth; Kalman Q = 0.005)

| Method | Energy removed |
|---|---:|
| Savitzky-Golay | 86 % |
| f-x decon | 64 % |
| SVD | 72 % |
| Kalman | 57 % |
| Kalman RTS | 64 % |
