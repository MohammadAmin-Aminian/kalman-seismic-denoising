function out = denoise_all(noisy, cfg, Q, dt, fx)
%DENOISE_ALL Apply every method; out{i} corresponds to METHOD_NAMES{i}.
%   Q, dt and the f-x settings default to the synthetic-gather configuration.
k = cfg.kalman;
if nargin < 3 || isempty(Q),  Q = k.Q; end
if nargin < 4 || isempty(dt), dt = cfg.synthetic.dt_s; end
if nargin < 5 || isempty(fx), fx = cfg.fx; end
out = cell(1, 5);
out{1} = savgol_smooth(noisy, cfg.savgol.order, cfg.savgol.frame);
out{2} = fx_decon_filter(noisy, dt, fx.operator_length, fx.mu_percent, fx.fmin_hz, fx.fmax_hz);
out{3} = svd_denoise(noisy, cfg.svd.rank);
out{4} = kalman1d(noisy, Q, k.R, k.P0, k.x0);
out{5} = kalman_rts(noisy, Q, k.R, k.P0, k.x0);
end
