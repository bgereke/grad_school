% Demo of interactive optimization of resolution enhancement for
% poorly resolved peaks.  Displays sliders for separate real-time control  
% of 2nd and 4th derivative weighting factors (factor and factor2) and 
% smooth width. Larger values of factor1 and factor2 will reduce the 
% peak width but will cause artifacts in the baseline near the peak. 
% Adjust the factors for the best trade-off for your application. 
% Use minimum smooth width needed to reduce excess noise. You can
% change the width of the peaks in line 24 and the peak shape in line 29.
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
noise=.01;
t=[1:8000];
amp=[1 2 3 4 5];  % Amplitudes of the peaks
pos=[300 400 500 600 700];   % Positions of the peaks
PeakWidth=162;
wid=[PeakWidth PeakWidth PeakWidth PeakWidth PeakWidth];   % Widths of the peaks
% A = matrix containing one unit-amplidude peak in each of its rows
A = zeros(length(pos),length(t));
for k=1:length(pos)
  A(k,:)=lorentzian(t,pos(k),wid(k)); % You can use any peak function here
end
% Multiply each row by the corresponding amplitude and adds them up
puresignal=amp*A;  
% Add noise
signal=pow(35,1:8000);%puresignal+noise.*randn(size(t));
SmoothWidth=round(PeakWidth/2.5);
factor1=round((PeakWidth.^2)./25);  % Second derivative weighting factor
factor2=round((PeakWidth.^4)./833);  % Fourth derivative weighting factor
% Plot the simulated signal
h=figure(1);
PlotRange=[SmoothWidth.*3:length(t)-SmoothWidth.*3];
h2=gca;
my=max(signal);
Enhancedsignal=enhance(signal,factor1,factor2,SmoothWidth);
Enhancedsignal(1:(SmoothWidth.*3))=signal(1:(SmoothWidth.*3));
Enhancedsignal(length(t)-SmoothWidth.*3:length(t))=signal(length(t)-SmoothWidth.*3:length(t));
plot(t,signal, t(PlotRange),Enhancedsignal(PlotRange),'r')
title(['factor1 = ' num2str(factor1)  '    factor2 = '  num2str(factor2) '    SmoothWidth= ' num2str(SmoothWidth)])
xlabel('BLUE = Original signal     RED = Resolution-enhanced signal')
my=max(signal);
axis([t(1) t(length(t)) -.1*my 2*my]);
grid on
% Add real-time sliders to graph for control of factor1 and factor2.
rtslid(h,@ResEnhanceF1,h2,1,'Scale',[0 4*factor1],'Def',factor1,'Back',[0.9 0.9 0.9],'Label','Factor 1','Position',[0.03 0.5 0.03 0.35]);
rtslid(h,@ResEnhanceF2,h2,0,'Scale',[0 4*factor2],'Def',factor2,'Back',[0.9 0.9 0.9],'Label','Factor 2','Position',[0.03 0.05 0.03 0.35]);
rtslid(h,@ResEnhanceF3,h2,0,'Scale',[0 200],'Def',SmoothWidth,'Back',[0.9 0.9 0.9],'Label','Smooth','Position',[0.95 0.1 0.03 0.35]);
