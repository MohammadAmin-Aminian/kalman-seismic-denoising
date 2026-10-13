function s = snr_db(clean, estimate)
%SNR_DB 10 log10( sum(clean^2) / sum((estimate - clean)^2) ).
e = estimate(:) - clean(:);
s = 10 * log10(sum(clean(:) .^ 2) / sum(e .^ 2));
end
