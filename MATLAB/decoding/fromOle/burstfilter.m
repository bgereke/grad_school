function [xb,xs] = burstfilter(x,dt);

xb = 0*x; 
xs = 0*x; 

j = 0;
k = 0;
for i=2:length(x)-1;
    if ((x(i) - x(i-1)) > dt) & ((x(i+1) - x(i)) > dt)   
        j = j + 1;
        xs(j) = x(i); 
    else
        k = k + 1;
        xb(k) = x(i); 
    end 
end
xs = xs(1:j);  
xb = xb(1:k);  

