function r = run_synthetic(cfg, sigma)
%RUN_SYNTHETIC Denoise the synthetic gather with noise level sigma.
%   r.snr and r.snr_orig are [input, METHOD_NAMES...] in dB.
[clean, t, h] = hyperbolic_gather(cfg.synthetic);
noisy = clean + sigma * load_noise(cfg.noise.file);
out = denoise_all(noisy, cfg);
snr = zeros(1, 6); snrOrig = zeros(1, 6);
snr(1) = snr_db(clean, noisy); snrOrig(1) = snr(1);
for i = 1:5
    snr(i + 1) = snr_db(clean, out{i});
    snrOrig(i + 1) = snr_db_original(clean, out{i});
end
r = struct('sigma', sigma, 't', t, 'h', h, 'clean', clean, 'noisy', noisy);
r.out = out;
r.snr = snr;
r.snr_orig = snrOrig;
end
