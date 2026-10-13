function h = gather_offsets(s)
%GATHER_OFFSETS Offsets (m) of the synthetic gather (row vector).
n = round((s.offset_stop_m - s.offset_start_m) / s.offset_step_m) + 1;
h = s.offset_start_m + s.offset_step_m * (0:n - 1);
if isfield(s, 'mirror') && s.mirror
    h = [-fliplr(h), h];
end
end
