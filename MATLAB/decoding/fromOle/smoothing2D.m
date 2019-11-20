function y = smooting2D(x,hx,hy)
% function y = smooting2D(x,hx,hy)
%

dx = ceil(5*hx);
dy=  ceil(5*hy);

for i=-dx:dx
    for j=-dy:dy   
        cx(dx+1+i,dy+1+j) = exp(-i^2/(2*hx^2) -j^2/(2*hy^2));
    end
end
cx = cx/sum(sum(cx));

y = conv2(x,cx,'same');

