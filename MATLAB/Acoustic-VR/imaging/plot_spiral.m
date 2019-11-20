function plot_spiral(trans,time,phase,id)

%produces spiral plot as seen in Danielson et al. 2016 Fig. 1C
%inputs:
%trans - binary vector of detected transients
%time - vector specifying time of each sample
%phase - vector specifying phase of each sample
%id - inner diameter of spiral (proportion of outer diamter)
%outputs: a figure

offset = max(time)/(1/id-1);
time = time + offset;
time = time/max(time);

plot(time.*cos(phase),time.*sin(phase),'-k','LineWidth',0.25);hold on
plot(time(trans==1).*cos(phase(trans==1)),time(trans==1).*sin(phase(trans==1)),'.r','MarkerSize',15)
axis square
axis off
hold off