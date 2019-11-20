function w = occupancy(p,tmin,tmax)
% w = occupancy(p,tmin,tmax)
% 
% How much time is spend at a certain location.


fprintf(1,'start: function occupancy, tmin=%d tmax=%d\n',tmin,tmax); 

bins = 150;
lgd  = 205; 
dt =   0.05; 
w = zeros(bins,1);  
j = 0;

ind = find(p(:,1) >tmin &  p(:,1) < tmax);
p = p(ind,:); 

for i=1:length(p)-1
    t = p(i,1) ;
    pos      = getpos(p,t);
    x        = p(i,2);
    y        = p(i,3);
    l        = getdist([x y]); 
    b        = round(l*bins/lgd); 
    if b > 0 & b <= bins
        w(b) = w(b) + dt;  
    end
end 

fprintf(1,'stop:  function occupancy\n'); 
