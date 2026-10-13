function B = savgol_coefficients(order, frame)
%SAVGOL_COEFFICIENTS Projection matrix of a least-squares polynomial fit.
if mod(frame, 2) == 0 || frame <= order
    error('savgol:badFrame', 'frame must be odd and larger than order');
end
m = (frame - 1) / 2;
s = (-m:m)' / m;                 % scaled positions for conditioning
V = s .^ (0:order);
B = V * ((V' * V) \ V');
end
