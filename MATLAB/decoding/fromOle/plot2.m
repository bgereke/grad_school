function plot2(x,S);

[C,l] = correlogram(x,x,0.2,0.001);
[Pxx,f] = psdd(C,1000);
subplot(2,1,1);
plot(l,C);
title(S) 
xlabel('sec')
ylabel('correlation')
subplot(2,1,2);
plot(f,log(Pxx));
xlabel('frequency (Hz)')
ylabel('power')
set(gca,'XLim',[0 140])

