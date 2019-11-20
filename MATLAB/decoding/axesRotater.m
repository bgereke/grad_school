function [xx, yy, theta] = axesRotater (x,y,varargin)

%varargin is for the angle. if it is already known and supplied

if isempty(varargin)
    length = abs(max(x) - min(x));
    cut = .2 * length;
    xb = x(x>min(x)+cut & x<max(x)-cut);
    yb = y(x>min(x)+cut & x<max(x)-cut);%try to only look at the data in the long part of the track to get the best line
    g = polyfit(xb,yb,1);
    %g(1)
    theta = atan(g(1));
else
    theta = varargin{1};
end

xx = x*cos(theta) + y*sin(theta);
yy = -x*sin(theta) + y*cos(theta);