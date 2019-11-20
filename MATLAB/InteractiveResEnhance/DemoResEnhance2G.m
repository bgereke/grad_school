% Demo of interactive optimization of derivative resolution enhancement
% for a single Gaussian peak. Displays sliders for separate real-time
% control of 2nd and 4th derivative weighting factors (factor and factor2)
% and smooth width. The peak width of the resolution-enhanced peak 
% is computed and displayed.  Larger values of factor1 and factor2 will
% reduce the peak width but will cause artifacts in the baseline near 
% the peak.  Adjust these two factors for the best compromise. 
% Use the minimum smooth width needed to reduce excess noise.
% Change the peak width in line 26 to see how the optimum factors change
% with peak width. Change the peak shape in line 28 to see how the
% optimum factors change with peak shape.
% Tom O'Haver, July 2006.  Slider function by Matthew Jones.

close
global t
global PlotRange
global signal
global Enhancedsignal
global SmoothWidth
global factor1
global factor2
format compact
clf;hold off
noise=.001;
t=[1:1000];
PeakWidth=200;
% Use any peak function here (gaussian, lorentzian, pearson, etc.)
puresignal = 10.*(gaussian(t,500,PeakWidth));
signal=puresignal+noise.*randn(size(t));
SmoothWidth=round(PeakWidth./10);
factor1=round((PeakWidth.^2)./25);  % Second derivative weighting factor
factor2=round((PeakWidth.^4)./833);  % Fourth derivative weighting factor
h=figure(1);
PlotRange=[SmoothWidth.*3:length(t)-SmoothWidth.*3];
plot(t(PlotRange),signal(PlotRange),'b')
h2=gca;
Enhancedsignal=enhance(signal,factor1,factor2,SmoothWidth);
y=Enhancedsignal';
x=t';
options = optimset('TolX',0.1);
start=[500 100];
estimated_lambda=FMINSEARCH('fitgauss',start,options,x,y,h);
MeasuredWidth= estimated_lambda(2);
PlotRange=[SmoothWidth.*3:length(t)-SmoothWidth.*3];
plot(t(PlotRange),signal(PlotRange), t(PlotRange),Enhancedsignal(PlotRange),'r')
title(['factor1 = ' num2str(factor1)  '    factor2 = '  num2str(factor2) '    SmoothWidth = ' num2str(SmoothWidth) '  PeakWidth = ' num2str(round(MeasuredWidth))])
xlabel('BLUE = Original signal     RED = Resolution-enhanced signal')
my=max(signal);
axis([t(1) t(length(t)) -.1*my 3*my]);
grid on

% Add real-time sliders to graph for control of factor1 and factor2.
rtslid(h,@ResEnhance2F1,h2,1,'Scale',[0 4*factor1],'Def',factor1,'Back',[0.9 0.9 0.9],'Label','Factor 1','Position',[0.03 0.5 0.03 0.35]);
rtslid(h,@ResEnhance2F2,h2,0,'Scale',[0 4*factor2],'Def',factor2,'Back',[0.9 0.9 0.9],'Label','Factor 2','Position',[0.03 0.05 0.03 0.35]);
rtslid(h,@ResEnhance2F3,h2,0,'Scale',[1 50],'Def',SmoothWidth,'Back',[0.9 0.9 0.9],'Label','Smooth','Position',[0.95 0.1 0.03 0.35]);

 