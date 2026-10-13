function xs = kalman_rts(z, Q, R, P0, x0)
%KALMAN_RTS Rauch-Tung-Striebel smoother for the KALMAN1D model.
%   Extension (not in the original study): forward filter followed by a
%   backward pass, which removes the phase delay of the causal filter.
if nargin < 3 || isempty(R),  R = 0.05; end
if nargin < 4 || isempty(P0), P0 = 1.01; end
if nargin < 5 || isempty(x0), x0 = 0;    end
rowIn = isrow(z);
if rowIn, z = z(:); end
[xf, ~, Pf] = kalman1d(z, Q, R, P0, x0);
xs = xf;
for k = size(xs, 1) - 1:-1:1
    Pp = Pf(k, :) + Q;                         % P_{k+1}^-
    C = zeros(size(Pp));
    pos = Pp > 0;
    C(pos) = Pf(k, pos) ./ Pp(pos);
    xs(k, :) = xf(k, :) + C .* (xs(k + 1, :) - xf(k, :));
end
if rowIn, xs = xs.'; end
end
