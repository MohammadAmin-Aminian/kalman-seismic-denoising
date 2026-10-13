function fd = run_field(cfg)
%RUN_FIELD Field record: AGC and denoising (no ground truth).
%   fd.removed is the fraction of energy removed by each method.
f = cfg.field; k = cfg.kalman;
shot = load_shot(f.file);
dt = shot.dt_s;
x = agc_gain(shot.d, agc_window_samples(f.agc_window_s, dt));
out = denoise_all(x, cfg, f.Q, dt, f.fx);
removed = zeros(1, 5);
for i = 1:5
    removed(i) = sum((x(:) - out{i}(:)) .^ 2) / sum(x(:) .^ 2);
end
ss = kalman_steady_state(f.Q, k.R);
fd = struct('t', (0:size(x, 1) - 1)' * dt, 'h', shot.offset_m, 'dt', dt, 'x', x, ...
            'removed', removed, 'K', ss.K, 'f3db_hz', cutoff_frequency(ss.K, dt), 'Q', f.Q);
fd.out = out;
end
