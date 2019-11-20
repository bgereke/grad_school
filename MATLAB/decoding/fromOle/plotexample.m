function plotexampple(xv,tv,pv,pfc,tspkc,tmin,dt)

close
subplot(2,1,1)
hold off

plot(xv,205*pv/150,'y-')
hold on 
plot(xv,205*pv/150,'rx')
plot(xv,205*tv/150,'k-'); 

set(gca,'XLim',[tmin tmin+dt]) 
set(gca,'YLim',[-5 210]) 
xlabel('Time (sec) ')
ylabel('Position (cm) ')
title('Bayesian Method, > 9 spikes')


subplot(2,1,2)
plot(xv,pfc,'k-');
xlabel('Time (sec) ')
ylabel('# of cells firing per dt ')
%hold on
%plot(xv,tspkc,'y-');
% plot(xv,tspkc./pfc,'k-');
set(gca,'XLim',[tmin tmin+dt]) 
