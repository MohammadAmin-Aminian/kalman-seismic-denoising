function [f, A] = amplitude_spectrum(d, dt)
%AMPLITUDE_SPECTRUM One-sided |FFT| averaged over traces (columns).
if isrow(d), d = d(:); end
nt = size(d, 1);
X = fft(double(d), [], 1);
nh = floor(nt / 2) + 1;
A = mean(abs(X(1:nh, :)), 2);
f = (0:nh - 1)' / (nt * dt);
end
