function p = smooting(x,h)
% function p = smooting(x,h)
%
% Perform Parzen smoothing using a Gaussian kernel of width h
% e.g. Bishop (1995), p. 52 
% NB. slow and can be optimized ... 

lgd = length(x); 
for i=1:lgd             

    sum = 0;
    for j=1:lgd              
        l1 = abs(i-j);
        l2 = lgd - abs(i-j);
        l = min(l1,l2); 
        if l < 5*h
            sum = sum + x(j)*exp(-l^2/(2*h^2))/(sqrt(2*pi)*h);      
        end
    end
%    p(i) = sum/lgd                 
    p(i) = sum;
end

