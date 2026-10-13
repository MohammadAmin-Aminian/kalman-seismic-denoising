function t = time_axis(dt, tmax)
%TIME_AXIS Column vector 0:dt:tmax (robust to rounding).
nt = floor(tmax / dt + 1e-9) + 1;
t = (0:nt - 1)' * dt;
end
