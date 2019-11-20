function [tv,xc]= interpsmooth(x,h)

dx = ceil(5*h);
i = -dx:dx;
cx = exp(-i.^2/(2*h^2) );
cx = cx/sum(sum(cx));


tv = min(x(2,:)):max(x(2,:));
xi = interp1(x(2,:),x(1,:),tv); 

xc = conv(xi,cx); 
xc = xc(dx:length(xc)-dx-1);
