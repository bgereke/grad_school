% Redraw the graph of DemoResEnhance when the any slider is changed
Enhancedsignal=enhance(signal,factor1,factor2,SmoothWidth);
Enhancedsignal(1:(SmoothWidth.*3))=signal(1:(SmoothWidth.*3));
Enhancedsignal(length(t)-SmoothWidth.*3:length(t))=signal(length(t)-SmoothWidth.*3:length(t));
plot(t,signal, t(PlotRange),Enhancedsignal(PlotRange),'r')
title(['factor1 = ' num2str(factor1)  '    factor2 = '  num2str(factor2) '    SmoothWidth= ' num2str(SmoothWidth)])
xlabel('BLUE = Original signal     RED = Resolution-enhanced signal')
my=max(signal);
axis([t(1) t(length(t)) -.1*my 2*my]);
grid on