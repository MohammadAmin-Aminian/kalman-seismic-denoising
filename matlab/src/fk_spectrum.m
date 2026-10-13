function [k, f, S] = fk_spectrum(d, dt, dx)
%FK_SPECTRUM Centred f-k amplitude spectrum in dB (k in 1/km, f >= 0).
[nt, nx] = size(d);
nft = 2 ^ nextpow2(nt);
nfx = 2 ^ nextpow2(nx);
D = fftshift(fft2(double(d), nft, nfx), 2);
D = abs(D(1:nft / 2 + 1, :));
S = 20 * log10(D / max(D(:)) + 1e-12);
f = (0:nft / 2)' / (nft * dt);
k = ((0:nfx - 1) - nfx / 2) / (nfx * dx) * 1000;
end
