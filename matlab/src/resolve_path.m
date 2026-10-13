function p = resolve_path(rel)
%RESOLVE_PATH Path relative to the repository root (absolute paths unchanged).
if isempty(rel) || rel(1) == '/' || rel(1) == '\' || (numel(rel) > 1 && rel(2) == ':')
    p = rel;
else
    p = fullfile(repo_root(), rel);
end
end
