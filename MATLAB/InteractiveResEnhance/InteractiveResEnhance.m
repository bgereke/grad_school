% Interactive optimization of derivative resolution enhancement for
% vector "signal". Resolution-enhanced signal is left in the vector 
% "Enhancedsignal". Displays sliders for separate real-time control  
% of 2nd and 4th derivative weighting factors (factor and factor2) and 
% smooth width. (Larger values of factor1 and factor2 will educe the 
% peak width but will cause artifacts in the baseline near 
% the peak.  Adjust the factors for the the best compromise.)
% Use the minimum smooth width needed to reduce excess noise. 
% If the range of the sliders is inappropriate for your signal,
% you can adjust the slider ranges in lines 27-29.
% Functions needed:ResEnhanceF1, ResEnhanceF2, ResEnhanceF3,
% rtslid, 
% Tom O'Haver, July 2006
% Slider function by Matthew Jones.
close
global t
global PlotRange
global signal
global SmoothWidth
global factor1
global factor2
format compact
clf;hold off
t=[1:length(signal)];
PeakWidth=length(signal)./30;
% Change the slider ranges in the next 3 lines
Factor1Range=2.*(PeakWidth.^2)/5;
Factor2Range=4.*(PeakWidth.^4)/370;
SmoothWidthRange=PeakWidth;
SmoothWidth=14;
factor1=0;
factor2=0;
% Plot the  signal
h=figure(1);
PlotRange=[SmoothWidth.*3:length(t)-SmoothWidth.*3];
plot(t(PlotRange),signal(PlotRange),'b')
h2=gca;
my=max(signal);
axis([t(1) t(length(t)) -.1*my 2*my]);
% Label the plot
title(['factor1 = ' num2str(factor1)  '    factor2 = '  num2str(factor2) '    SmoothWidth= ' num2str(SmoothWidth)])
xlabel('BLUE = Original signal     RED = Resolution-enhanced signal')
grid on
% Add real-time sliders to graph for control of factor1 and factor2.
rtslid(h,@ResEnhanceF1,h2,1,'Scale',[0 Factor1Range],'Def',0,'Back',[0.9 0.9 0.9],'Label','Factor 1','Position',[0.03 0.5 0.03 0.35]);
rtslid(h,@ResEnhanceF2,h2,0,'Scale',[0 Factor2Range],'Def',0,'Back',[0.9 0.9 0.9],'Label','Factor 2','Position',[0.03 0.05 0.03 0.35]);
rtslid(h,@ResEnhanceF3,h2,0,'Scale',[0 SmoothWidthRange],'Def',SmoothWidth,'Back',[0.9 0.9 0.9],'Label','Smooth','Position',[0.95 0.1 0.03 0.35]);
