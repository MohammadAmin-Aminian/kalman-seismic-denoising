function cfg = load_config(path)
%LOAD_CONFIG Read the shared JSON configuration (config/default.json).
if nargin < 1 || isempty(path)
    path = fullfile(repo_root(), 'config', 'default.json');
end
cfg = jsondecode(fileread(path));
end
