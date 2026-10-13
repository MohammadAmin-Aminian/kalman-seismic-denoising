function [d, t, h] = hyperbolic_gather(s)
%HYPERBOLIC_GATHER Clean gather of hyperbolic reflections.
%   d(t,h) = sum_j amp_j * ricker(t - sqrt(tau_j^2 + (h/v_j)^2)).
t = time_axis(s.dt_s, s.tmax_s);
h = gather_offsets(s);
d = zeros(numel(t), numel(h));
for j = 1:numel(s.tau_s)
    tt = sqrt(s.tau_s(j) ^ 2 + (h / s.v_mps(j)) .^ 2);
    d = d + s.amp(j) * ricker_wavelet(t - tt, s.f0_hz);
end
end
