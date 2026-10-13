function r = repo_root()
%REPO_ROOT Repository root folder (two levels above this file).
r = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
