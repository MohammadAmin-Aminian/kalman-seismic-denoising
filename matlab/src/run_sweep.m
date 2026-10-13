function sw = run_sweep(cfg)
%RUN_SWEEP Output SNR of the Kalman filter/smoother as a function of Q.
k = cfg.kalman;
Q = logspace(cfg.sweep.log10_Q_min, cfg.sweep.log10_Q_max, cfg.sweep.n);
clean = hyperbolic_gather(cfg.synthetic);
unit = load_noise(cfg.noise.file);
sig = cfg.noise.sigma_levels(:)';
sw = struct('Q', Q, 'sigma', sig);
sw.K = zeros(size(Q));
for i = 1:numel(Q)
    ss = kalman_steady_state(Q(i), k.R);
    sw.K(i) = ss.K;
end
sw.forward = zeros(numel(sig), numel(Q));
sw.rts = zeros(numel(sig), numel(Q));
sw.input = zeros(1, numel(sig));
for s = 1:numel(sig)
    noisy = clean + sig(s) * unit;
    sw.input(s) = snr_db(clean, noisy);
    for i = 1:numel(Q)
        sw.forward(s, i) = snr_db(clean, kalman1d(noisy, Q(i), k.R, k.P0, k.x0));
        sw.rts(s, i) = snr_db(clean, kalman_rts(noisy, Q(i), k.R, k.P0, k.x0));
    end
end
end
