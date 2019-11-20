function plotphasedist(i,xmin,xmax,ymax)

div = 5; 
load(strcat('nall',num2str(div)));
subplot(3,1,1)
plotsplit(nallp,i,div,ymax); 
set(gca,'XLim',[xmin xmax]) 
legend('1','2','3','4','5')

div = 7; 
load(strcat('nall',num2str(div)));
subplot(3,1,2)
plotsplit(nallp,i,div,ymax); 
set(gca,'XLim',[xmin xmax]) 
legend('1','2','3','4','5','6','7')

div = 9; 
load(strcat('nall',num2str(div)));
subplot(3,1,3)
plotsplit(nallp,i,div,ymax); 
set(gca,'XLim',[xmin xmax]) 
legend('1','2','3','4','5','6','7','8','9')
