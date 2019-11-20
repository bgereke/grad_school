function ResEnhanceF3(n,h)
% Redraw the graph when the SmoothWidth slider is changed
% Tom O'Haver, July 2006
global t
global PlotRange
global signal
global Enhancedsignal
global SmoothWidth
global factor1
global factor2
SmoothWidth=1+round(n);
PlotRange=[SmoothWidth.*3:length(t)-SmoothWidth.*3];
ResEnhanceRedraw
  
