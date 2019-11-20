function [xb,xs] = burstfilter(x,dt);

xb = 0*x; 
xs = 0*x; 

j = 0;
k = 0;
for i=2:length(x);
    if (x(i) - x(i-1)) < dt   
        j = j + 1;
        xb(j) = x(i); 
    else
        k = k + 1;
        xs(k) = x(i); 
    end 
end
xb = xb(1:j);  
xs = xs(1:k);  

