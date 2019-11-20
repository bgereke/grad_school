function ResEnhance2F3(n,h)
% Re-draw graph when SmoothWidth slider is changed
% Tom O'Haver, July 2006
global t
global PlotRange
global signal
global Enhancedsignal
global SmoothWidth
global factor1
global factor2
SmoothWidth=1+round(n);
ResEnhance2Redraw