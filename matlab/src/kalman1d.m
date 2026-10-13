function [x, K, P] = kalman1d(z, Q, R, P0, x0)
%KALMAN1D Forward (causal) random-walk Kalman filter along dimension 1.
%   [x, K, P] = KALMAN1D(z, Q, R, P0, x0) filters every column of z
%   independently with the scalar model
%       x_k = x_{k-1} + w_k,  w_k ~ N(0, Q)
%       z_k = x_k + v_k,      v_k ~ N(0, R)
%   Defaults: R = 0.05, P0 = 1.01, x0 = 0. NaN samples are treated as
%   missing observations (prediction only). Returns the filtered estimate
%   x, the Kalman gain K and the posterior variance P (same size as z).
if nargin < 3 || isempty(R),  R = 0.05; end
if nargin < 4 || isempty(P0), P0 = 1.01; end
if nargin < 5 || isempty(x0), x0 = 0;    end
if ~(Q >= 0 && R > 0 && P0 >= 0)
    error('kalman1d:badParameters', 'require Q >= 0, R > 0 and P0 >= 0');
end
rowIn = isrow(z);
if rowIn, z = z(:); end
z = double(z);
[nt, ntr] = size(z);
x = zeros(nt, ntr); K = zeros(nt, ntr); P = zeros(nt, ntr);
xk = repmat(x0, 1, ntr);
Pk = repmat(P0, 1, ntr);
for k = 1:nt
    if k > 1
        Pm = Pk + Q;
    else
        Pm = Pk;            % the prior already describes the first sample
    end
    obs = ~isnan(z(k, :));
    Kk = zeros(1, ntr);
    Kk(obs) = Pm(obs) ./ (Pm(obs) + R);
    innov = zeros(1, ntr);
    innov(obs) = z(k, obs) - xk(obs);
    xk = xk + Kk .* innov;
    Pk = (1 - Kk) .* Pm;
    x(k, :) = xk; K(k, :) = Kk; P(k, :) = Pk;
end
if rowIn
    x = x.'; K = K.'; P = P.';
end
end
