function plotsplit(x,i,div,ymax)

cells = 38; 
hold on
for j=1:div 
        k = j + div *(i - 1); 

        plot((1:150)*205/150,x(:,1+(i-1)*div : div*i )); 
end

line(205*[25 25]/150,[0 ymax])
line(205*[75 75]/150,[0 ymax])
line(205*[125 125]/150,[0 ymax])

xlabel('position (bins)')
ylabel('firing rate (1/s)')
