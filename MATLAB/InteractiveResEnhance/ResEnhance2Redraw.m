% Re-draw graph for DemoResEnhance2 when when any slider is changed
axes(h);
Enhancedsignal=enhance(signal,factor1,factor2,SmoothWidth);
y=Enhancedsignal';
x=t';
options = optimset('TolX',0.1);
start=[500 50];
estimated_lambda=FMINSEARCH('fitgauss',start,options,x,y,h);
MeasuredWidth= estimated_lambda(2);
PlotRange=[SmoothWidth.*3:length(t)-SmoothWidth.*3];
plot(t(PlotRange),signal(PlotRange), t(PlotRange),Enhancedsignal(PlotRange),'r')
title(['factor1 = ' num2str(factor1)  '    factor2 = '  num2str(factor2) '    SmoothWidth = ' num2str(SmoothWidth) '  PeakWidth = ' num2str(round(MeasuredWidth))])
xlabel('BLUE = Original peak     RED = Resolution-enhanced peak')
my=max(signal);
axis([t(1) t(length(t)) -.1*my 3*my]);
grid on