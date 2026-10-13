function ch = run_chirp(cfg)
%RUN_CHIRP Measure the filter gain on a linear chirp and compare with theory.
c = cfg.chirp; k = cfg.kalman;
[x, t, finst] = chirp_signal(c.dt_s, c.duration_s, c.f_start_hz, c.f_end_hz);
noise = load_noise(c.noise_file);
noisy = x + c.sigma * noise(:, 1);
yClean = kalman1d(x, k.Q, k.R, k.P0, k.x0);
yNoisy = kalman1d(noisy, k.Q, k.R, k.P0, k.x0);
ss = kalman_steady_state(k.Q, k.R);
gainMeas = abs(analytic_signal(yClean)) ./ abs(analytic_signal(x));
gainTheory = abs(ewma_response(ss.K, finst, c.dt_s));
n = numel(t);
mid = (floor(0.1 * n) + 1):floor(0.9 * n);     % same samples as Python slice
ch = struct('t', t, 'f_inst', finst, 'clean', x, 'noisy', noisy, ...
            'filtered_clean', yClean, 'filtered_noisy', yNoisy, ...
            'gain_measured', gainMeas, 'gain_theory', gainTheory, ...
            'max_gain_error', max(abs(gainMeas(mid) - gainTheory(mid))), ...
            'K', ss.K, 'f3db_hz', cutoff_frequency(ss.K, c.dt_s));
end
