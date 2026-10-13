function s = load_shot(rel)
%LOAD_SHOT Shot gather d(nt, ntraces), offsets (m), sampling interval (s).
S = load(resolve_path(rel));
s = struct('d', double(S.d), 'offset_m', double(S.offset_m(:))', ...
           'dt_s', double(S.dt_s), 'source', strtrim(char(S.source)));
end
