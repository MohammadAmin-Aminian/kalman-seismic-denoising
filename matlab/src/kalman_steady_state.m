function ss = kalman_steady_state(Q, R)
%KALMAN_STEADY_STATE Steady-state variances and gain of KALMAN1D.
%   P^- = (Q + sqrt(Q^2 + 4 Q R)) / 2,  K = P^- / (P^- + R),  P^+ = (1-K) P^-.
if nargin < 2 || isempty(R), R = 0.05; end
if ~(Q >= 0 && R > 0)
    error('kalman_steady_state:badParameters', 'require Q >= 0 and R > 0');
end
Pm = 0.5 * (Q + sqrt(Q * Q + 4 * Q * R));
K = Pm / (Pm + R);
ss = struct('P_prior', Pm, 'P_post', (1 - K) * Pm, 'K', K);
end
