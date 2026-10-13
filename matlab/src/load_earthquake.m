function eq = load_earthquake(rel)
%LOAD_EARTHQUAKE Vertical broadband record (counts) and metadata.
S = load(resolve_path(rel));
eq = struct('counts', double(S.counts(:)), 'fs_hz', double(S.fs_hz), ...
            'station', strtrim(char(S.station)), 'start_time', strtrim(char(S.start_time)));
end
