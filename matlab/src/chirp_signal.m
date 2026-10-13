function [x, t, finst] = chirp_signal(dt, duration, fStart, fEnd)
%CHIRP_SIGNAL Unit-amplitude linear chirp (column vectors).
t = time_axis(dt, duration);
rate = (fEnd - fStart) / duration;
x = sin(2 * pi * (fStart * t + 0.5 * rate * t .^ 2));
finst = fStart + rate * t;
end
