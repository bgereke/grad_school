function x = getratio(p,t)
% function x = getratio(p,t)
% p = [time thetaper ratio]                  
% x = [ratio Ttheta]             

if t < p(1,1)
    error('Error in getpos.m. t smaller than smallest t in p ...')  
end
if t > p(length(p),1)
    error('Error in getpos.m. t bigger than largest in p ...')  
end

i = floor((t - p(1,1))/0.05)-100;  
if i < 1
    i = 1;
end
if i > length(p)
    i = length(p); 
end

while p(i,1) > t
   i = i - 50; 
end

while p(i,1) <=  t
    i = i+1;
end
i=i-1; 
r1 = p(i,3);
r2 = p(i+1,3);

x = [r1 + (r2 - r1)*(t - p(i,1))/(p(i+1,1) - p(i,1)) p(i+1,2)] ; 
