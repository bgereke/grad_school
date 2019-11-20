function ta = acf(x,y,p)

ta = zeros(p+1,1) ;
N = max(size(y)) ;
y = y-mean(y);
x = x-mean(x);
yvar = y'*y;
xvar = x'*x;
den = sqrt(xvar*yvar);

for i = 0:p   
   ta(i+1) = x(1:N-i)'*y(i+1:N)/den;
end







