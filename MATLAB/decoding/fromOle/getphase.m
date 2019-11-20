function x = getphase(theta,t)
% function x = getphase(theta,t)
% theta : [times of theta peak]
% 
% Find the phase of firing at time t using the information about
% theta             

Tmax = 0.167;
Tmin = 0.100;

i  = locate(theta,t); 

T = theta(i+1) - theta(i); 
x = (t - theta(i))/(theta(i+1) - theta(i));
if T > Tmax | T < Tmin
    x = -1; 
end
