function n = load_noise(rel)
%LOAD_NOISE Unit-variance Gaussian noise stored in single precision.
S = load(resolve_path(rel), 'noise');
n = double(S.noise);
end
