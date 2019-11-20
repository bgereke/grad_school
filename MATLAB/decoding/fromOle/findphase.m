function d = findphase(p)

global x;

sum = 0; 
x0 = min(x(:,1)); 
for i = 1:length(x)
    xp = x(i,1);
    yp = x(i,2);

    yl = p(1)*(xp-x0) + p(2);

    d1 = mod(abs(yp - yl),1);      
    d2 = mod(abs(yp - (yl+1)),1);      
    d3 = mod(abs(yp - (yl-1)),1);      

    d = min([d1 d2 d3]);
    sum = sum+d*d; 
end
d = sum; 
