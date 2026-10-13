function H = ewma_response(K, f, dt)
%EWMA_RESPONSE Frequency response of the steady-state filter.
%   H(f) = K / (1 - (1 - K) exp(-i 2 pi f dt)).
w = 2 * pi * f * dt;
H = K ./ (1 - (1 - K) .* exp(-1i * w));
end
