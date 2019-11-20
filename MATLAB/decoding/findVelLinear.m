function [vel] = findVelLinear(xx,t)
%rotate separately

vel = zeros(1,length(xx));

for i = 2:1:(length(xx)-1)
    vel(i) = (xx(i+1)-xx(i-1))/(t(i+1)-t(i-1)); 
end
vel(length(xx)) = (xx(length(xx))-xx(length(xx)-1))/(t(length(t))-t(length(t)-1));



