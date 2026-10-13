function z = analytic_signal(x)
%ANALYTIC_SIGNAL Analytic signal via FFT (equivalent to HILBERT, base MATLAB).
x = double(x(:));
n = numel(x);
h = zeros(n, 1);
h(1) = 1;
if mod(n, 2) == 0
    h(n / 2 + 1) = 1;
    h(2:n / 2) = 2;
else
    h(2:(n + 1) / 2) = 2;
end
z = ifft(fft(x) .* h);
end
