function d = getdist(l)
% function d = getdist(l)
%
% l = [x1 x2],  location
% 
 
% define the corners
c1 = [38 192];
c2 = [132 66];
c3 = [195 209];
 
% define lengths of arms
d1 = norm(c2 - c1);
d2 = norm(c3 - c2);
d3 = norm(c1 - c3);

a = getarm(l);

if a == 1
    v1 = l  - c1;
    v2 = c2 - c1;     
    vp = v2*(v1*v2')/(norm(v2)^2);  
    d  = norm(vp);  
end
if a == 2

    v1 = l  - c2;
    v2 = c3 - c2;     
    vp = v2*(v1*v2')/(norm(v2)^2);  
    d  = norm(vp) + d1;  
end
if a == 3
    v1 = l  - c3;
    v2 = c1 - c3;     
    vp = v2*(v1*v2')/(norm(v2)^2);  
    d  = norm(vp) + d1 + d2;  
end

% 0.4348 cm/pixel
d = 0.4348*d;  
