function w = ricker_wavelet(t, f0)
%RICKER_WAVELET Ricker wavelet with peak frequency f0 (Hz), unit peak at t = 0.
a = (pi * f0 * t) .^ 2;
w = (1 - 2 * a) .* exp(-a);
end
