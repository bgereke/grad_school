function plotres(z,xmin,xmax);


plot(z(:,1),z(:,2),'.')
hold on
plot(z(:,1),z(:,2)+1,'.')
hold off
set(gca,'XLim',[xmin xmax]);
set(gca,'YLim',[0 2]);
title('phase'); 
ylabel('phase'); 
xlabel('position (cm)'); 

orient tall
