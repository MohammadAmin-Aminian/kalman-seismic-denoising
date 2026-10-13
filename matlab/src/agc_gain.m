function y = agc_gain(d, window)
%AGC_GAIN Automatic gain control: divide by the RMS in a centred window.
%   y = AGC_GAIN(d, window) works on every column of d. The window (odd
%   number of samples) shrinks at the trace ends. The gather is first
%   scaled by its maximum absolute value (AGC is scale invariant); a
%   stabilisation of 1e-6 times the largest RMS avoids division by zero.
if mod(window, 2) == 0 || window < 1
    error('agc_gain:badWindow', 'window must be a positive odd number of samples');
end
rowIn = isrow(d);
if rowIn, d = d(:); end
d = double(d);
peak = max(abs(d(:)));
if peak == 0
    y = d;
    if rowIn, y = y.'; end
    return
end
X = d / peak;
nt = size(X, 1);
h = (window - 1) / 2;
c = [zeros(1, size(X, 2)); cumsum(X .^ 2, 1)];
i = (1:nt)';
lo = max(i - h, 1);
hi = min(i + h, nt);
ms = (c(hi + 1, :) - c(lo, :)) ./ (hi - lo + 1);
rms = sqrt(max(ms, 0));
y = X ./ (rms + 1e-6 * max(rms(:)));
if rowIn, y = y.'; end
end
