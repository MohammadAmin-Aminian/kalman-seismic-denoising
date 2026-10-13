classdef TestKalmanDenoise < matlab.unittest.TestCase
    %TESTKALMANDENOISE Unit, known-answer and Python-parity tests (base MATLAB).
    %   Run from the repository root:
    %       results = runtests('matlab/tests');

    properties
        Ref
        Cfg
    end

    methods (TestClassSetup)
        function setPath(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(here, '..', 'src')));
            tc.Ref = load(fullfile(here, '..', '..', 'tests', 'reference', 'python_reference.mat'));
            tc.Cfg = load_config();
        end
    end

    methods (Test)
        % ---------------------------------------------------- Kalman filter
        function gainConvergesToClosedForm(tc)
            for Q = [1e-5, 4e-3, 0.5]
                [~, K] = kalman1d(zeros(4000, 1), Q, 0.05);
                ss = kalman_steady_state(Q, 0.05);
                tc.verifyLessThan(abs(K(end) - ss.K), 1e-12);
            end
        end

        function constantIsRecovered(tc)
            x = kalman1d(3.7 * ones(500, 1), 1e-3, 0.05);
            tc.verifyLessThan(abs(x(end) - 3.7), 1e-6);
        end

        function steadyStateIsFirstOrderRecursion(tc)
            z = sin((1:3000)' / 7) + cos((1:3000)' / 3);
            Q = 4e-3; R = 0.05;
            x = kalman1d(z, Q, R);
            ss = kalman_steady_state(Q, R);
            y = x(1000);
            for k = 1001:3000
                y = (1 - ss.K) * y + ss.K * z(k);
            end
            tc.verifyLessThan(abs(y - x(end)), 1e-9);
        end

        function mseMatchesPosteriorVariance(tc)
            s = RandStream('mt19937ar', 'Seed', 1);
            Q = 0.01; R = 0.05; n = 4000; m = 200;
            truth = cumsum(sqrt(Q) * randn(s, n, m), 1);
            z = truth + sqrt(R) * randn(s, n, m);
            x = kalman1d(z, Q, R, 1e6);
            mse = mean((x(501:end, :) - truth(501:end, :)) .^ 2, 'all');
            ss = kalman_steady_state(Q, R);
            tc.verifyLessThan(abs(mse / ss.P_post - 1), 0.03);
            xs = kalman_rts(z, Q, R, 1e6);
            mseS = mean((xs(501:end - 500, :) - truth(501:end - 500, :)) .^ 2, 'all');
            tc.verifyLessThan(mseS, 0.8 * mse);
        end

        function missingSamplesAreSkipped(tc)
            z = sin(linspace(0, 10, 400))';
            z(101:110) = NaN;
            [x, K] = kalman1d(z, 4e-3);
            tc.verifyTrue(all(isfinite(x)));
            tc.verifyEqual(K(101:110), zeros(10, 1));
            tc.verifyEqual(x(101:110), repmat(x(100), 10, 1));
        end

        function rowAndColumnInputs(tc)
            z = sin((1:200) / 5);
            tc.verifyEqual(kalman1d(z, 4e-3), kalman1d(z(:), 4e-3).');
        end

        function cutoffIsMinus3dB(tc)
            ss = kalman_steady_state(4e-3, 0.05);
            fc = cutoff_frequency(ss.K, 0.002);
            tc.verifyLessThan(abs(abs(ewma_response(ss.K, fc, 0.002)) - 1 / sqrt(2)), 1e-12);
            tc.verifyLessThan(abs(fc - 22.5836), 1e-3);
        end

        function invalidParametersError(tc)
            tc.verifyError(@() kalman1d(zeros(10, 1), -1, 0.05), 'kalman1d:badParameters');
            tc.verifyError(@() kalman1d(zeros(10, 1), 1e-3, 0), 'kalman1d:badParameters');
        end

        % ---------------------------------------------------- baselines
        function savgolReproducesPolynomials(tc)
            t = linspace(-1, 1, 200)';
            p = 1 - 2 * t + 0.5 * t .^ 3 - t .^ 4;
            tc.verifyLessThan(max(abs(savgol_smooth(p, 4, 27) - p)), 1e-10);
        end

        function fxImprovesPlaneWaves(tc)
            dt = 0.002; t = (0:1000)' * dt; x = (0:59) * 10;
            d = ricker_wavelet(t - (0.3 + 4e-4 * x), 25) ...
                - 0.8 * ricker_wavelet(t - (0.9 - 2e-4 * x), 25) ...
                + 0.6 * ricker_wavelet(t - 1.4 - 0 * x, 25);
            s = RandStream('mt19937ar', 'Seed', 0);
            noisy = d + 0.2 * randn(s, size(d));
            out = fx_decon_filter(noisy, dt, 8, 1, 0, 100);
            tc.verifyTrue(isreal(out));
            tc.verifyGreaterThan(snr_db(d, out) - snr_db(d, noisy), 8);
        end

        function svdRecoversRankOne(tc)
            A = (1:50)' * cos(0:19);
            tc.verifyLessThan(max(abs(svd_denoise(A, 1) - A), [], 'all'), 1e-10);
        end

        function snrDefinitions(tc)
            c = sin(linspace(0, 20, 1000));
            tc.verifyLessThan(abs(snr_db(c, 0.5 * c) - 20 * log10(2)), 1e-12);
            tc.verifyLessThan(abs(snr_db_original(c, 0.5 * c)), 1e-12);
        end

        % ---------------------------------------------------- parity with Python
        function parityKalmanOnChirp(tc)
            k = tc.Cfg.kalman;
            [x, K, P] = kalman1d(tc.Ref.chirp_noisy(:), k.Q, k.R, k.P0, k.x0);
            tc.verifyLessThan(max(abs(x - tc.Ref.chirp_kf(:))), 1e-12);
            tc.verifyLessThan(max(abs(K - tc.Ref.chirp_K(:))), 1e-14);
            tc.verifyLessThan(max(abs(P - tc.Ref.chirp_P(:))), 1e-14);
            xs = kalman_rts(tc.Ref.chirp_noisy(:), k.Q, k.R, k.P0, k.x0);
            tc.verifyLessThan(max(abs(xs - tc.Ref.chirp_rts(:))), 1e-12);
        end

        function paritySyntheticGatherAndMethods(tc)
            sig = tc.Cfg.noise.sigma_levels;
            r = run_synthetic(tc.Cfg, sig(end));
            cols = double(tc.Ref.cols);
            tc.verifyLessThan(max(abs(r.clean(:, cols) - tc.Ref.clean_cols), [], 'all'), 1e-12);
            tc.verifyLessThan(max(abs(r.noisy(:, cols) - tc.Ref.noisy_cols), [], 'all'), 1e-12);
            for i = 1:5
                ref = tc.Ref.(sprintf('out%d_cols', i));
                err = max(abs(r.out{i}(:, cols) - ref), [], 'all') / max(abs(ref), [], 'all');
                tc.verifyLessThan(err, 1e-8, sprintf('method %d', i));
            end
            tc.verifyLessThan(max(abs(r.snr - tc.Ref.snr_hi(:)')), 1e-6);
            tc.verifyLessThan(max(abs(r.snr_orig - tc.Ref.snr_orig_hi(:)')), 1e-6);
            r2 = run_synthetic(tc.Cfg, sig(1));
            tc.verifyLessThan(max(abs(r2.snr - tc.Ref.snr_lo(:)')), 1e-6);
        end

        function paritySpectra(tc)
            [d, ~, ~] = hyperbolic_gather(tc.Cfg.synthetic);
            [f, A] = amplitude_spectrum(d, tc.Cfg.synthetic.dt_s);
            tc.verifyLessThan(max(abs(f - tc.Ref.amp_f(:))), 1e-9);
            tc.verifyLessThan(max(abs(A - tc.Ref.amp_A(:))) / max(A), 1e-10);
            [k, ff, S] = fk_spectrum(d, tc.Cfg.synthetic.dt_s, tc.Cfg.synthetic.offset_step_m);
            tc.verifyLessThan(max(abs(k(:) - tc.Ref.fk_k(:))), 1e-9);
            tc.verifyLessThan(max(abs(ff(:) - tc.Ref.fk_f(:))), 1e-9);
            S0 = S(:, numel(k) / 2 + 1);
            % compare above -100 dB only: below that level the values are set by
            % floating-point rounding, which differs between FFT libraries
            R0 = tc.Ref.fk_S_col_k0(:);              % column, same shape as S0
            ok = R0 > -100;
            tc.verifyLessThan(max(abs(S0(ok) - R0(ok))), 1e-6);
        end

        function paritySweep(tc)
            sw = run_sweep(tc.Cfg);
            tc.verifyLessThan(max(abs(sw.Q - tc.Ref.sweep_Q(:)') ./ sw.Q), 1e-12);
            tc.verifyLessThan(max(abs(sw.forward(1, :) - tc.Ref.sweep_fwd_lo(:)')), 1e-6);
            tc.verifyLessThan(max(abs(sw.forward(2, :) - tc.Ref.sweep_fwd_hi(:)')), 1e-6);
            tc.verifyLessThan(max(abs(sw.rts(1, :) - tc.Ref.sweep_rts_lo(:)')), 1e-6);
            tc.verifyLessThan(max(abs(sw.rts(2, :) - tc.Ref.sweep_rts_hi(:)')), 1e-6);
        end

        function parityChirpAndEarthquake(tc)
            ch = run_chirp(tc.Cfg);
            tc.verifyLessThan(abs(ch.max_gain_error - tc.Ref.chirp_gain_error), 1e-9);
            tc.verifyLessThan(ch.max_gain_error, 0.05);
            eq = run_earthquake(tc.Cfg);
            tc.verifyLessThan(max(abs([eq.main.K, eq.alt.K] - tc.Ref.eq_K(:)')), 1e-14);
            tc.verifyLessThan(max(abs([eq.main.rms_ratio, eq.alt.rms_ratio] - tc.Ref.eq_rms_ratio(:)')), 1e-10);
            tc.verifyLessThan(max(abs([eq.main.f3db_hz, eq.alt.f3db_hz] - tc.Ref.eq_f3db(:)')), 1e-10);
            yh = eq.main.y(1:2000);
            tc.verifyLessThan(max(abs(yh - tc.Ref.eq_y_main_head(:))) / max(abs(yh)), 1e-10);
        end

        function reproduceAllRuns(tc)
            out = tempname;
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(here, '..')));
            res = evalc('r = reproduce_all(''figures'', false, ''outdir'', out);'); %#ok<NASGU>
            tc.verifyTrue(isfile(fullfile(out, 'results', 'metrics_matlab.json')));
            tc.verifyTrue(isfile(fullfile(out, 'results', 'snr_table_matlab.md')));
        end
    end
end
