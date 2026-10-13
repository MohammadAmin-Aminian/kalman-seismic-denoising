function mo = run_modelled(cfg)
%RUN_MODELLED Modelled shot gather: AGC, add Gaussian noise, denoise, score.
%   mo.snr is [input, METHOD_NAMES...] in dB against the noise-free gather.
m = cfg.modelled; k = cfg.kalman;
shot = load_shot(m.file);
dt = shot.dt_s;
clean = agc_gain(shot.d, agc_window_samples(m.agc_window_s, dt));
rmsv = sqrt(mean(clean(:) .^ 2));
noisy = clean + m.sigma * rmsv * load_noise(m.noise_file);
out = denoise_all(noisy, cfg, m.Q, dt, m.fx);
snr = zeros(1, 6);
snr(1) = snr_db(clean, noisy);
for i = 1:5
    snr(i + 1) = snr_db(clean, out{i});
end
Qs = logspace(cfg.sweep.log10_Q_min, cfg.sweep.log10_Q_max, cfg.sweep.n);
fwd = zeros(size(Qs)); rts = zeros(size(Qs));
for i = 1:numel(Qs)
    fwd(i) = snr_db(clean, kalman1d(noisy, Qs(i), k.R, k.P0, k.x0));
    rts(i) = snr_db(clean, kalman_rts(noisy, Qs(i), k.R, k.P0, k.x0));
end
mo = struct('t', (0:size(clean, 1) - 1)' * dt, 'h', shot.offset_m, 'dt', dt, ...
            'clean', clean, 'noisy', noisy, 'snr', snr, 'sweep_Q', Qs, ...
            'sweep_forward', fwd, 'sweep_rts', rts, 'Q', m.Q, 'sigma', m.sigma);
mo.out = out;
end
