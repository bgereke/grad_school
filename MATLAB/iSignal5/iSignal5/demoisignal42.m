% Self-running demo script of iSignal 4.2
% Calls plotit verison 6 
x=0:.005:2;y=humps(x);Data=[x;y];
disp('The built-in Matlab Humps function, displayed in iSignal')
isignal(Data,0.6,4,0,0,1,0,0,0,0,0,0);
subplot(2,1,1)
title('The built-in Matlab Humps function, displayed in iSignal')
pause(2)
isignal(Data,0.245,0.615,4,3,1,0,0,0,0,0,0,0);
subplot(2,1,1)
title('Viewing the first peak near x=0.4')
drawnow
pause(2)
isignal(Data,0.945,0.735,4,3,1,0,0,0,0,0,0);
subplot(2,1,1)
title('Viewing the second peak near x=0.9')
pause(2)
isignal(Data,0.9,0.5,1,3,1,4);
subplot(2,1,1)
title('4th derivative of the peak at x=0.9')
pause(2)
isignal(Data,0.3,0.5,1,3,1,0,1,220,5400);
subplot(2,1,1)
title('Peak sharpening applied to the peak at x=0.3')
pause(2)
clf
x=-5:.01:5;
y=exp(-(x).^2);
for n=1:1000,
    if randn()>2,y(n)=rand()+y(n);
    end
end
isignal([x;y],0,20,0,3,0,0,0,10,1000,0,0,0);
subplot(2,1,1)
title('Smooth signal with many sharp spikes')
pause(2)
isignal([x;y],0,20,0,3,0,0,0,10,1000,0,1,0);
subplot(2,1,1)
title('Spikes removed using median filter (M key)')
pause(2)
x=0:.1:60; y=sin(x)+sin(10.*x);
[pY,PowerSpectrum]=isignal([x;y],30,30,4,3,1,0,0,1,0,0,0,1);
pause(2)
clf
plot(PowerSpectrum(:,1),PowerSpectrum(:,2))
title('Power spectrum plotted on linear coordinates')
pause(2)
loglog(PowerSpectrum(:,1),PowerSpectrum(:,2))
title('Power spectrum plotted on log-log coordinates')
pause(2)
% Demonstration script for iSignal function. It generates a test signal
% consisting of 3 peaks on a curved baseline, then runs iSignal
%   T. C. O'Haver, December 2013, January 2015
increment=.1;
x=1:increment:700;
% For each simulated peak, compute the amplitude, position, and width
% Peak 1 is the baseline
amp=[100 .5 2 3 3.5];  % Amplitudes of the peaks  (CHANGE if desired)
pos=[-200 100 250 400 600];   % Positions of the peaks (CHANGE if desired)
wid=[700 50 50 50 50];   % Widths of the peaks (CHANGE if desired)
Noise=.01; % Amount of random noise added to the signal. (CHANGE if desired) 
% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
      % Create a series of peaks of different x-positions
      A(k,:)=gaussian(x,pos(k),wid(k)); % Gaussian or Lorentzian peaks
      % A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      p=p+1;
end
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));  % Optionally adds random noise
DataMatrix=[x' y']; % Assembles x and y vectors into data matrix
disp('-----------------------------------------------------------------')
disp('Detection and measurement of 4 peaks heights 1, 2, 3, 3.5, with ')
disp('equal widths, superimposed in a very strong curved baseline.')
disp('The objective is to extract a measure that is proportional to ')
disp('the peak height but independent of the baseline strength. In this demo, ')
disp('differentiation (with smoothing) is used to supress the baseline.')
isignal(DataMatrix);
disp('Signal consisting of 4 weak peaks on a strong curved baseline.')
subplot(211);
title('Signal consisting of 4 weak peaks on a strong curved baseline.')
pause(2)
isignal(DataMatrix,250.5,499,0,1,0,0,0,10,1000,0,0,0);
disp('Press Ctrl-A to zoom out to see entire signal')
subplot(211);
title('Press Ctrl-A to zoom out to see entire signal')
pause(2)
isignal(DataMatrix,250.5,499,0,1,0,1,0,10,1000,0,0,0);
disp('Press D key to compute first derivative')
subplot(211);
title('Press D key to compute first derivative')
pause(2)
isignal(DataMatrix,250.5,499,0,1,0,2,0,10,1000,0,0,0);
disp('Press D key again to compute second derivative')
subplot(211);
title('Press D key again to compute second derivative')
pause(2)
isignal(DataMatrix,250.5,499,1,12,0,2,0,10,1000,0,0,0);
subplot(211);
title('Second derivative with 12 point rectangular smooth')
pause(1)
isignal(DataMatrix,250.5,499,1,25,0,2,0,10,1000,0,0,0);
subplot(211);
title('Second derivative with 25 point rectangular smooth')
pause(1)
isignal(DataMatrix,250.5,499,1,51,0,2,0,10,1000,0,0,0);
subplot(211);
title('Second derivative with 51 point rectangular smooth')
pause(1)
isignal(DataMatrix,250.5,499,1,95,0,2,0,10,1000,0,0,0);
subplot(211);
title('Second derivative with 95 point rectangular smooth')
pause(1)
isignal(DataMatrix,250.5,499,2,95,0,2,0,10,1000,0,0,0);
subplot(211);
title('Second derivative with 95 point triangular smooth')
pause(1)
isignal(DataMatrix,250.5,499,3,95,0,2,0,10,1000,0,0,0);
subplot(211);
title('Second derivative with 95 point Gaussian smooth')
pause(1)
isignal(DataMatrix,250.5,499,4,95,0,2,0,10,1000,0,0,0);
subplot(211);
title('Second derivative with 95 point Savitsky-Golay smooth')
pause(1)
isignal(DataMatrix,250.5,499,3,95,0,2,0,10,1000,0,0,0);
subplot(211);
title('Second derivative with 95 point Gaussian smooth')
pause(2)
[pY,PowerSpectrum,maxy,miny,area,stdev]=isignal(DataMatrix,100,147.4,3,95,0,2,0,10,1000,0,0,0);
s1=maxy-miny;
subplot(211);
title('Measure peak-to-peak signal of peak 1')
pause(2)
[pY,PowerSpectrum,maxy,miny,area,stdev]=isignal(DataMatrix,250,147.4,3,95,0,2,0,10,1000,0,0,0);
s2=maxy-miny;
subplot(211);
title('Measure peak-to-peak signal of peak 2')
pause(2)
[pY,PowerSpectrum,maxy,miny,area,stdev]=isignal(DataMatrix,400,147.4,3,95,0,2,0,10,1000,0,0,0);
s3=maxy-miny;
subplot(211);
title('Measure peak-to-peak signal of peak 3')
pause(2)
[pY,PowerSpectrum,maxy,miny,area,stdev]=isignal(DataMatrix,600,147.4,3,95,0,2,0,10,1000,0,0,0);
s4=maxy-miny;
subplot(211);
title('Measure peak-to-peak signal of peak 4')
pause(2)
clf
% Note: following line requires plotit version 6 (January 2015)
% Download from http://terpconnect.umd.edu/~toh/spectrum/plotit.m
plotit([amp(2) amp(3) amp(4) amp(5)],[s1 s2 s3 s4],1,'or');    
title('Calibration plot of derivative signal ampliude against true peak amplitude')  
xlabel('True amplitude of the 4 original peaks')
ylabel('Peak-to-peak derivative amplitude of each peak')

disp('End of demo.')
