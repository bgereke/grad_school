function Enhancedsignal=enhance(signal,factor1,factor2,SmoothWidth)
% Resolution enhancement function by derivative method. the
% arguments factor1 and factor 2 are 2nd and 4th derivative weighting
% factors. Larger values of factor1 and factor2 will reduce the 
% peak width but will cause artifacts in the baseline near 
% the peak.  Adjust the factors for the the best compromise. 
% Use minimum smooth width needed to reduce excess noise. 
% Use InteracticeResEnhance.m to adjust these factors 
% interactively on your own signals.
% See DemoResEnhance.m and DemoResEnhance2.m for examples of use.
% Functions required: secderiv.m, tsmooth.m
d2=secderiv(signal);  % Computes second derivative
d4=secderiv(d2);   % Computes fourth derivative
Enhancedsignal = signal-factor1.*tsmooth(d2,SmoothWidth)+...
factor2.*tsmooth(tsmooth(tsmooth(d4,SmoothWidth),SmoothWidth),SmoothWidth);