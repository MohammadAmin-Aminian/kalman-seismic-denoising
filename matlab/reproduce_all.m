function results = reproduce_all(varargin)
%REPRODUCE_ALL Regenerate all numbers and figures of the study (base MATLAB).
%   REPRODUCE_ALL() writes results/metrics_matlab.json,
%   results/snr_table_matlab.md and figures/matlab/*.png.
%
%   REPRODUCE_ALL('figures', false)      numbers only
%   REPRODUCE_ALL('outdir', folder)      write outputs below another folder
%   REPRODUCE_ALL('config', file)        use another JSON configuration
%
%   Run from the repository root or from the matlab/ folder:
%       >> cd matlab; reproduce_all

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'src'));

p = inputParser;
p.addParameter('figures', true);
p.addParameter('outdir', '');
p.addParameter('config', '');
p.parse(varargin{:});
opt = p.Results;

cfg = load_config(opt.config);
root = opt.outdir;
if isempty(root), root = repo_root(); end
resDir = fullfile(root, 'results');
figDir = fullfile(root, 'figures', 'matlab');
if ~exist(resDir, 'dir'), mkdir(resDir); end

sig = cfg.noise.sigma_levels(:)';
names = method_names();
syn = cell(1, numel(sig));
for s = 1:numel(sig)
    syn{s} = run_synthetic(cfg, sig(s));
end
sw = run_sweep(cfg);
ch = run_chirp(cfg);
eq = run_earthquake(cfg);

% ---------------------------------------------------------------- numbers
snrTable = zeros(numel(sig), 6);
snrOrig = zeros(numel(sig), 6);
for s = 1:numel(sig)
    snrTable(s, :) = syn{s}.snr;
    snrOrig(s, :) = syn{s}.snr_orig;
end
[bestF, iF] = max(sw.forward, [], 2);
[bestR, iR] = max(sw.rts, [], 2);
results = struct();
results.columns = [{'input'}, names];
results.sigma = sig;
results.synthetic_snr_db = snrTable;
results.synthetic_snr_db_original_definition = snrOrig;
results.sweep = struct('Q', sw.Q, 'forward', sw.forward, 'rts', sw.rts, ...
    'best_forward_Q', sw.Q(iF), 'best_forward_snr_db', bestF', ...
    'best_rts_Q', sw.Q(iR), 'best_rts_snr_db', bestR');
results.chirp = struct('K', ch.K, 'f3db_hz', ch.f3db_hz, 'max_gain_error', ch.max_gain_error);
results.earthquake = struct( ...
    'main', rmfield(eq.main, 'y'), 'alt', rmfield(eq.alt, 'y'));

fid = fopen(fullfile(resDir, 'metrics_matlab.json'), 'w');
fprintf(fid, '%s\n', jsonencode(results));
fclose(fid);

rows = [{'Input (noisy)'}, names, {'Kalman, best Q (sweep)', 'Kalman RTS, best Q (sweep)'}];
vals = [snrTable, bestF, bestR];                       % sigma x rows
fid = fopen(fullfile(resDir, 'snr_table_matlab.md'), 'w');
fprintf(fid, 'Output SNR (dB), synthetic gather (MATLAB).\n\n| Method |');
fprintf(fid, ' sigma = %g |', sig);
fprintf(fid, '\n|---|%s\n', repmat('---:|', 1, numel(sig)));
for i = 1:numel(rows)
    fprintf(fid, '| %s |', rows{i});
    fprintf(fid, ' %.2f |', vals(:, i));
    fprintf(fid, '\n');
end
fclose(fid);
type(fullfile(resDir, 'snr_table_matlab.md'));

if ~opt.figures, return; end
if ~exist(figDir, 'dir'), mkdir(figDir); end

% ---------------------------------------------------------------- figures
C = struct('ink', [11 11 11] / 255, 'muted', [138 137 134] / 255, ...
    'kal', hex('2a78d6'), 'sg', hex('eb6834'), 'fx', hex('1baf7a'), ...
    'svd', hex('eda100'), 'rts', hex('6250d6'));
mcol = {C.sg, C.fx, C.svd, C.kal, C.rts};
lo = syn{1}; hi = syn{end};
dt = cfg.synthetic.dt_s;
clip = 0.5 * max(abs(lo.clean(:)));
gray_r = flipud(gray(256));

% fig01 gathers
f = newfig(9, 4.2);
panels = {lo.clean, '(a) Clean'; lo.noisy, sprintf('(b) Noisy, \\sigma = %g (SNR %.1f dB)', lo.sigma, lo.snr(1)); ...
          hi.noisy, sprintf('(c) Noisy, \\sigma = %g (SNR %.1f dB)', hi.sigma, hi.snr(1))};
for i = 1:3
    subplot(1, 3, i);
    gatherImage(panels{i, 1}, lo.t, lo.h, clip, gray_r);
    title(panels{i, 2});
    if i == 1, ylabel('Time (s)'); end
end
save_png(f, fullfile(figDir, 'fig01_synthetic_gathers.png'));

% fig02/03 methods
pairs = {hi, 'fig02_methods_sigma_high.png'; lo, 'fig03_methods_sigma_low.png'};
for q = 1:2
    r = pairs{q, 1};
    f = newfig(11, 6.4);
    for i = 1:5
        subplot(2, 5, i);
        gatherImage(r.out{i}, r.t, r.h, clip, gray_r);
        title({names{i}, sprintf('SNR %.1f dB', r.snr(i + 1))});
        if i == 1, ylabel('Time (s)'); end
        subplot(2, 5, 5 + i);
        gatherImage(r.noisy - r.out{i}, r.t, r.h, clip, gray_r);
        title('removed (noisy - output)');
        if i == 1, ylabel('Time (s)'); end
    end
    save_png(f, fullfile(figDir, pairs{q, 2}));
end

% fig04 trace comparison
j = cfg.display.trace_index;
sel = hi.t >= 0.35 & hi.t <= 1.25;
f = newfig(7.5, 6.2);
tr = {hi.noisy(:, j), 'Noisy', C.muted; hi.out{4}(:, j), 'Kalman', C.kal; ...
      hi.out{5}(:, j), 'Kalman RTS', C.rts; hi.out{1}(:, j), 'Savitzky-Golay', C.sg};
for i = 1:4
    subplot(4, 1, i);
    plot(hi.t(sel), hi.clean(sel, j), 'Color', C.ink, 'LineWidth', 1); hold on
    plot(hi.t(sel), tr{i, 1}(sel), 'Color', tr{i, 3}, 'LineWidth', 1.2);
    legend({'Clean', tr{i, 2}}, 'Location', 'northeast', 'Box', 'off');
    ylabel('Amplitude'); ylim([-2.5 4.2]); box off
    if i == 1, title(sprintf('Trace %d (offset %.0f m), \\sigma = %g', j, hi.h(j), hi.sigma)); end
end
xlabel('Time (s)');
save_png(f, fullfile(figDir, 'fig04_trace_comparison.png'));

% fig05 amplitude spectra
f = newfig(6.5, 3.8);
[fr, A] = amplitude_spectrum(hi.clean, dt); semilogy(fr, A, 'Color', C.ink, 'LineWidth', 1.6); hold on
[fr, A] = amplitude_spectrum(hi.noisy, dt); semilogy(fr, A, 'Color', C.muted);
order = [4 5 1 2];
for i = order
    [fr, A] = amplitude_spectrum(hi.out{i}, dt); semilogy(fr, A, 'Color', mcol{i}, 'LineWidth', 1.1);
end
xlim([0 cfg.display.fmax_hz]); ylim([1e-2 5e2]); box off
xlabel('Frequency (Hz)'); ylabel('Mean amplitude spectrum');
title(sprintf('Amplitude spectra, \\sigma = %g', hi.sigma));
legend([{'Clean', 'Noisy'}, names(order)], 'NumColumns', 3, 'Box', 'off');
save_png(f, fullfile(figDir, 'fig05_amplitude_spectra.png'));

% fig06 f-k spectra
f = newfig(9, 3.6);
fkp = {hi.clean, '(a) Clean'; hi.noisy, '(b) Noisy'; hi.out{4}, '(c) Kalman'};
for i = 1:3
    subplot(1, 3, i);
    [kk, ff, S] = fk_spectrum(fkp{i, 1}, dt, cfg.synthetic.offset_step_m);
    imagesc(kk, ff, S, [-40 0]); axis xy; ylim([0 cfg.display.fmax_hz]);
    colormap(gca, blues());
    title(fkp{i, 2}); xlabel('Wavenumber (1/km)');
    if i == 1, ylabel('Frequency (Hz)'); end
    if i == 3, cb = colorbar; cb.Label.String = 'Amplitude (dB)'; end
end
save_png(f, fullfile(figDir, 'fig06_fk_spectra.png'));

% fig07 SNR vs Q
f = newfig(9, 3.5);
for s = 1:numel(sig)
    subplot(1, numel(sig), s);
    semilogx(sw.Q, sw.forward(s, :), '-o', 'Color', C.kal, 'MarkerSize', 3, 'MarkerFaceColor', C.kal); hold on
    semilogx(sw.Q, sw.rts(s, :), '-s', 'Color', C.rts, 'MarkerSize', 3, 'MarkerFaceColor', C.rts);
    semilogx(sw.Q([1 end]), sw.input(s) * [1 1], '--', 'Color', C.muted);
    yl = ylim; semilogx(cfg.kalman.Q * [1 1], yl, '-', 'Color', [.76 .76 .72]); ylim(yl);
    xlabel('Process-noise variance Q (R = 0.05)'); ylabel('Output SNR (dB)');
    title(sprintf('\\sigma = %g', sig(s))); box off
    if s == 1
        legend({'Kalman (forward)', 'Kalman RTS', 'Input (no filtering)', 'Q = 0.004 (original)'}, ...
               'Location', 'west', 'Box', 'off');
    end
end
save_png(f, fullfile(figDir, 'fig07_snr_vs_Q.png'));

% fig08 filter response
f = newfig(7.5, 3.6);
fr = linspace(0, 0.5 / dt, 1000);
nt = 4096;
w = ricker_wavelet(((0:nt - 1) - nt / 2) * dt, cfg.synthetic.f0_hz);
W = abs(fft(w)); W = W(1:nt / 2 + 1); fw = (0:nt / 2) / (nt * dt);
area(fw, W / max(W), 'FaceColor', [.91 .91 .89], 'EdgeColor', 'none'); hold on
Qs = [1e-4, cfg.kalman.Q, 0.05, 0.5];
pal = {C.kal, C.sg, C.fx, C.rts};
labels = {sprintf('Ricker %g Hz (normalised)', cfg.synthetic.f0_hz)};
for i = 1:4
    ss = kalman_steady_state(Qs(i), cfg.kalman.R);
    plot(fr, abs(ewma_response(ss.K, fr, dt)), 'Color', pal{i}, 'LineWidth', 1.2);
    labels{end + 1} = sprintf('Q = %g (K = %.3f)', Qs(i), ss.K); %#ok<AGROW>
end
xlim([0 0.5 / dt]); ylim([0 1.05]); box off
xlabel('Frequency (Hz)'); ylabel('|H(f)| of steady-state filter');
title('The random-walk Kalman filter is a first-order low-pass filter');
legend(labels, 'Location', 'eastoutside', 'Box', 'off');
save_png(f, fullfile(figDir, 'fig08_filter_response.png'));

% fig09 chirp validation
f = newfig(7.5, 5.8);
win = [0 0.6; 3.4 3.6];
for i = 1:2
    subplot(2, 2, i);
    s2 = ch.t >= win(i, 1) & ch.t <= win(i, 2);
    plot(ch.t(s2), ch.noisy(s2), 'Color', C.muted, 'LineWidth', 0.6); hold on
    plot(ch.t(s2), ch.clean(s2), 'Color', C.ink, 'LineWidth', 1.0);
    plot(ch.t(s2), ch.filtered_noisy(s2), 'Color', C.kal, 'LineWidth', 1.1);
    ylim([-2.4 2.4]); xlabel('Time (s)'); box off
    fa = interp1(ch.t, ch.f_inst, win(i, :));
    title(sprintf('(%s) %g-%g s (%.0f-%.0f Hz)', 'a' + i - 1, win(i, 1), win(i, 2), fa(1), fa(2)));
    if i == 1
        ylabel('Amplitude');
        legend({'Chirp + noise', 'Chirp', 'Kalman output'}, 'Location', 'southwest', 'Box', 'off');
    end
end
subplot(2, 2, [3 4]);
s3 = ch.f_inst >= 5 & ch.f_inst <= 115;
plot(ch.f_inst(s3), ch.gain_measured(s3), 'Color', C.kal, 'LineWidth', 1.4); hold on
plot(ch.f_inst(s3), ch.gain_theory(s3), '--', 'Color', C.ink);
plot(ch.f3db_hz * [1 1], [0 1.05], '-', 'Color', [.76 .76 .72]);
ylim([0 1.05]); box off
xlabel('Instantaneous frequency (Hz)'); ylabel('Amplitude gain');
title(sprintf('(c) Gain on the noise-free chirp (K = %.3f, -3 dB at %.1f Hz)', ch.K, ch.f3db_hz));
legend({'Measured (envelope ratio)', 'Theory |K/(1-(1-K)e^{-i\omega\Deltat})|'}, 'Box', 'off');
save_png(f, fullfile(figDir, 'fig09_chirp_validation.png'));

% fig10 earthquake
f = newfig(8, 6.4);
subplot(3, 1, 1);
plot(eq.t, eq.x, 'Color', C.muted, 'LineWidth', 0.5); hold on
plot(eq.t, eq.alt.y, 'Color', C.sg, 'LineWidth', 0.6);
legend({'Raw (demeaned)', sprintf('Kalman Q = %g', eq.alt.Q)}, 'Box', 'off', 'NumColumns', 2);
ylabel('Counts'); xlim([eq.t(1) eq.t(end)]); box off
title(sprintf('%s vertical record (40 Hz)', eq.meta.station));
subplot(3, 1, 2);
plot(eq.t, eq.x, 'Color', C.muted, 'LineWidth', 0.5); hold on
plot(eq.t, eq.main.y, 'Color', C.kal, 'LineWidth', 0.7);
legend({'Raw (demeaned)', sprintf('Kalman Q = %g', eq.main.Q)}, 'Box', 'off', 'NumColumns', 2);
ylabel('Counts'); xlim([eq.t(1) eq.t(end)]); box off
xlabel(sprintf('Time since %s UTC (s)', eq.meta.start_time(1:19)));
subplot(3, 1, 3);
sigs = {eq.x, C.muted, 'Raw'; eq.alt.y, C.sg, sprintf('Q = %g', eq.alt.Q); eq.main.y, C.kal, sprintf('Q = %g', eq.main.Q)};
for i = 1:3
    [fr, A] = amplitude_spectrum(sigs{i, 1}, eq.dt);
    loglog(fr(2:end), A(2:end), 'Color', sigs{i, 2}, 'LineWidth', 0.8); hold on
end
yl = ylim;
loglog(eq.alt.f3db_hz * [1 1], yl, '--', 'Color', C.sg);
loglog(eq.main.f3db_hz * [1 1], yl, '--', 'Color', C.kal);
xlim([1e-2 0.5 / eq.dt]); box off
xlabel('Frequency (Hz)  (dashed: -3 dB cut-off of each filter)'); ylabel('Amplitude');
legend(sigs(:, 3), 'Location', 'southwest', 'Box', 'off');
save_png(f, fullfile(figDir, 'fig10_earthquake.png'));
fprintf('figures written to %s\n', figDir);
end

% ------------------------------------------------------------- helpers
function c = hex(s)
c = [hex2dec(s(1:2)), hex2dec(s(3:4)), hex2dec(s(5:6))] / 255;
end

function f = newfig(wIn, hIn)
f = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [1 1 wIn hIn], ...
           'PaperPositionMode', 'auto');
set(f, 'DefaultAxesFontSize', 9, 'DefaultAxesTickDir', 'out', 'DefaultLineLineWidth', 1.2);
end

function gatherImage(d, t, h, clip, cmap)
imagesc(h / 1000, t, d, [-clip clip]);
colormap(gca, cmap);
xlabel('Offset (km)');
box off
end

function m = blues()
a = [247 251 255; 198 219 239; 107 174 214; 33 113 181; 8 48 107] / 255;
x = linspace(0, 1, size(a, 1));
xi = linspace(0, 1, 256)';
m = [interp1(x, a(:, 1), xi), interp1(x, a(:, 2), xi), interp1(x, a(:, 3), xi)];
end

function save_png(f, path)
print(f, path, '-dpng', '-r200');
close(f);
end
