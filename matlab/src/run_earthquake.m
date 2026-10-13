function eq = run_earthquake(cfg)
%RUN_EARTHQUAKE Apply the filter to the IU.DAV vertical record (demeaned counts).
e = cfg.earthquake;
m = load_earthquake(e.file);
x = m.counts - mean(m.counts);
dt = 1 / m.fs_hz;
eq = struct('meta', m, 'dt', dt, 'x', x, 't', (0:numel(x) - 1)' * dt);
names = {'main', 'alt'};
Qs = [e.Q, e.Q_alt];
for i = 1:2
    y = kalman1d(x, Qs(i), e.R, cfg.kalman.P0, 0);
    ss = kalman_steady_state(Qs(i), e.R);
    eq.(names{i}) = struct('Q', Qs(i), 'y', y, 'K', ss.K, ...
        'f3db_hz', cutoff_frequency(ss.K, dt), 'rms_ratio', sqrt(mean(y .^ 2) / mean(x .^ 2)));
end
end
