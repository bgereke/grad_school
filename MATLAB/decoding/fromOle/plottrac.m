function plottrac(p)
% function plottrac(p)
% 
% Uses p to plot the track of the rat using the function
% getpos.m to find location. 

x = zeros(20000,2); 
i = 1;
for t=1e4:0.05:1.1e4   
    x(i,:) = getpos(p,t); 
    i = i + 1;
end
plot(x(:,1),x(:,2),'.')

