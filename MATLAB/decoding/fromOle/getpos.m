function x = getpos(p,t)
% function x = getpos(p,t)
% p = [time posx posy angle]
% x = [x1 x2 angle]
% Estimate the position for time t based on the data in p

if t < p(1,1)
     % error('Error in getpos.m. t smaller than smallest t in p ...')  
     x = [-1 -1 -1];
     return
end
if t > p(length(p),1)
     % error('Error in getpos.m. t bigger than largest in p ...')  
     x = [-1 -1 -1];
     return
end


i = locate(p,t);

x1 = p(i,2);
x2 = p(i+1,2);
y1 = p(i,3);
y2 = p(i+1,3);
z1 = p(i,4);
z2 = p(i+1,4);

x(1) = x1 + (x2 - x1)*(t - p(i,1))/(p(i+1,1) - p(i,1)); 
x(2) = y1 + (y2 - y1)*(t - p(i,1))/(p(i+1,1) - p(i,1));
x(3) = z1 + (z2 - z1)*(t - p(i,1))/(p(i+1,1) - p(i,1));
