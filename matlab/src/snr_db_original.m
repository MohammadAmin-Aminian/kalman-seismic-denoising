function s = snr_db_original(clean, estimate)
%SNR_DB_ORIGINAL SNR as computed in the original scripts (for comparison only).
%   10 log10( sum(estimate^2) / sum((estimate - clean)^2) ) - flawed: the
%   numerator uses the estimate, not the clean signal.
e = estimate(:) - clean(:);
s = 10 * log10(sum(estimate(:) .^ 2) / sum(e .^ 2));
end
