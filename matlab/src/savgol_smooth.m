function y = savgol_smooth(x, order, frame)
%SAVGOL_SMOOTH Savitzky-Golay smoothing along dimension 1 (base MATLAB).
%   Same edge rule as SGOLAYFILT: the first/last (frame-1)/2 samples are
%   evaluated from the polynomial fitted to the first/last frame.
rowIn = isrow(x);
if rowIn, x = x(:); end
x = double(x);
[nt, ntr] = size(x);
if nt < frame, error('savgol:short', 'signal shorter than the frame'); end
B = savgol_coefficients(order, frame);
m = (frame - 1) / 2;
c = B(m + 1, :);
y = zeros(nt, ntr);
for i = 1:frame
    y(m + 1:nt - m, :) = y(m + 1:nt - m, :) + c(i) * x(i:nt - frame + i, :);
end
y(1:m, :) = B(1:m, :) * x(1:frame, :);
y(nt - m + 1:nt, :) = B(m + 2:end, :) * x(nt - frame + 1:nt, :);
if rowIn, y = y.'; end
end
