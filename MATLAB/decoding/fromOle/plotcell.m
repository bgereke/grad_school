function [x1,x2] = plotcell(x)

tmin = 1e4;
tmax = 1.1e4;      

load pos

x1 = size(length(x),1);
x2 = size(length(x),1);

for i=1:length(x);
    w =getpos(p,x(i));
    x1(i) = w(1);
    x2(i) = w(2);
%    plot(x1(i),x2(i),'.');
%    if mod(i,10) == 0
%         drawnow;
%    end

end


 
