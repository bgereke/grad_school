function DemoResEnhanceFactor2(n,h)
% Re-draw graph when factor2 slider is changed
% Tom O'Haver, July 2006
global t
global PlotRange
global signal
global Enhancedsignal
global SmoothWidth
global factor1
global factor2
factor2=round(n);
ResEnhance2Redraw
