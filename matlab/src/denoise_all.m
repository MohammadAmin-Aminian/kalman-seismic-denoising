function out = denoise_all(noisy, cfg, Q)
%DENOISE_ALL Apply every method; out{i} corresponds to METHOD_NAMES{i}.
k = cfg.kalman;
if nargin < 3 || isempty(Q), Q = k.Q; end
fx = cfg.fx;
out = cell(1, 5);
out{1} = savgol_smooth(noisy, cfg.savgol.order, cfg.savgol.frame);
out{2} = fx_decon_filter(noisy, cfg.synthetic.dt_s, fx.operator_length, fx.mu_percent, ...
                         fx.fmin_hz, fx.fmax_hz);
out{3} = svd_denoise(noisy, cfg.svd.rank);
out{4} = kalman1d(noisy, Q, k.R, k.P0, k.x0);
out{5} = kalman_rts(noisy, Q, k.R, k.P0, k.x0);
end
