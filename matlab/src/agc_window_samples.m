function n = agc_window_samples(window_s, dt)
%AGC_WINDOW_SAMPLES Odd window length (samples) closest to window_s.
n = 2 * round(window_s / dt / 2) + 1;
end
