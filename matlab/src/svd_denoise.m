function out = svd_denoise(d, r)
%SVD_DENOISE Rank-r reconstruction (sum of the first r eigenimages).
[U, S, V] = svd(double(d), 'econ');
out = U(:, 1:r) * S(1:r, 1:r) * V(:, 1:r)';
end
