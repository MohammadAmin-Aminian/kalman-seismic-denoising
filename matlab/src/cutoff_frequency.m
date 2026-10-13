function fc = cutoff_frequency(K, dt)
%CUTOFF_FREQUENCY -3 dB frequency (Hz) of the steady-state filter (NaN if none).
a = 1 - K;
if a <= 0
    fc = NaN; return
end
c = (1 + a * a - 2 * K * K) / (2 * a);
if c < -1
    fc = NaN; return
end
fc = acos(min(c, 1)) / (2 * pi * dt);
end
