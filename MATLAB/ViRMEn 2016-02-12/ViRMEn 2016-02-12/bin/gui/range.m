function r = range(a,d)
% range function, for users who don't have the statistics toolbox

if nargin == 1
    r = max(a) - min(a);
else
    r = max(a,[],d) - min(a,[],d);
end
