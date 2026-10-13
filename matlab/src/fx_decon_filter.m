function out = fx_decon_filter(d, dt, lf, mu, fmin, fmax)
%FX_DECON_FILTER f-x deconvolution (Canales, 1984) - base MATLAB.
%   For each frequency in [fmin, fmax] forward and backward spatial
%   prediction filters of length lf are estimated by damped least squares
%   (damping mu percent of the zero-lag autocorrelation); predictions are
%   averaged. Frequencies outside [fmin, fmax] are set to zero.
if nargin < 3 || isempty(lf),   lf = 8;   end
if nargin < 4 || isempty(mu),   mu = 0;   end
if nargin < 5 || isempty(fmin), fmin = 0; end
if nargin < 6 || isempty(fmax), fmax = 0.5 / dt; end
d = double(d);
[nt, ntr] = size(d);
nf = 2 ^ nextpow2(nt);
ilow = max(floor(fmin * dt * nf), 0);
ihigh = min(floor(fmax * dt * nf), nf / 2);
D = fft(d, nf, 1);
F = complex(zeros(nf, ntr));
Bk = complex(zeros(nf, ntr));
for k = ilow:ihigh
    [yf, yb] = ar_predict(D(k + 1, :).', lf, mu);
    F(k + 1, :) = yf.';
    Bk(k + 1, :) = yb.';
end
for k = nf / 2 + 1:nf - 1          % Hermitian symmetry -> real output
    F(k + 1, :) = conj(F(nf - k + 1, :));
    Bk(k + 1, :) = conj(Bk(nf - k + 1, :));
end
of = real(ifft(F, [], 1));
ob = real(ifft(Bk, [], 1));
out = of(1:nt, :) + ob(1:nt, :);
out(:, lf + 1:ntr - lf) = out(:, lf + 1:ntr - lf) / 2;
end

function [yf, yb] = ar_predict(x, lf, mu)
n = numel(x);
yf = complex(zeros(n, 1));
yb = complex(zeros(n, 1));
if n <= lf, return; end
Mf = complex(zeros(n - lf, lf));
Mb = complex(zeros(n - lf, lf));
for j = 1:lf
    Mf(:, j) = x(lf - j + 1:n - j);     % x(i-j) predicts x(i)
    Mb(:, j) = x(j + 1:n - lf + j);     % x(i+j) predicts x(i)
end
A = Mf' * Mf;
beta = real(A(1, 1)) * mu / 100;
a = (A + beta * eye(lf)) \ (Mf' * x(lf + 1:n));
yf(lf + 1:n) = Mf * a;
A = Mb' * Mb;
beta = real(A(1, 1)) * mu / 100;
b = (A + beta * eye(lf)) \ (Mb' * x(1:n - lf));
yb(1:n - lf) = Mb * b;
end
