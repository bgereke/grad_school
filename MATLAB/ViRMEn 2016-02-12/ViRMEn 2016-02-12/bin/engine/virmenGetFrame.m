function M = virmenGetFrame(w)
% M = virmenGetFrame(w)
%   Obtains the current image displayed in ViRMEn window w.
%   Output is a height x width x 3 matrix of RGB values

M = virmenOpenGLRoutines(5,w);
M = permute(M,[3 2 1]);