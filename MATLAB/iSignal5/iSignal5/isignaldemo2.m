function isignaldemo2
% Demonstration script for iSignal function. It generates a test signal
% consisting of 3 peaks on a curved baseline, then runs iSignal
%   T. C. O'Haver, December 2013
increment=.1;
x=[1:increment:500];
pausetime=.1;

% For each simulated peak, compute the amplitude, position, and width
% Peak 1 is the baseline
amp=[100 1 2 3];  % Amplitudes of the peaks  (Change if desired)
pos=[-200 100 250 400];   % Positions of the peaks (Change if desired)
wid=[700 50 50 50];   % Widths of the peaks (Change if desired)
Noise=.01; % Amount of random noise added to the signal. (Change if desired) 

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
demodata=[x' y']; % Assembles x and y vectors into data matrix


disp('-----------------------------------------------------------------')
disp('Detection and measurement of three peaks heights 1, 2, and 3, but with ')
disp('equal widths, superimposed in a very strong curved baseline.')
disp('The objective is to extract a measure that is proportional to ')
disp('the peak height but independent of the baseline strength')
disp('Suggested approaches: (a) Use baseline subtraction to remove the baseline;')
disp('or (b) use differentiation (with smoothing) to supress the baseline.')

isignal(demodata);
disp('Press Ctrl-A to zoom out to see entire signal')
subplot(211);title('Press Ctrl-A to zoom out to see entire signal')
pause(pausetime)
isignal(demodata,250.5,499,0,1,0,0,0,10,1000,0,0,0);
disp('Press Ctrl-A to zoom out to see entire signal')
subplot(211);title('Press Ctrl-A to zoom out to see entire signal')
pause(pausetime)
isignal(demodata,250.5,499,0,1,0,1,0,10,1000,0,0,0);
disp('Press D key to compute first derivative')
subplot(211);title('Press D key to compute first derivative')
pause(pausetime)
isignal(demodata,250.5,499,0,1,0,2,0,10,1000,0,0,0);
disp('Press D key again to compute second derivative')
subplot(211);title('Press D key again to compute second derivative')
pause(pausetime)
disp('End of demo.')

% ----------------------------------------------------------------------
function g = gaussian(x,pos,wid)
%  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
%  X may be scalar, vector, or matrix, pos and wid both scalar
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)).^2);

% ----------------------------------------------------------------------
function g = lorentzian(x,position,width)
% lorentzian(x,position,width) Lorentzian function.
% where x may be scalar, vector, or matrix
% position and width scalar
% T. C. O'Haver, 1988
% Example: lorentzian([1 2 3],2,2) gives result [0.5 1 0.5]
g=ones(size(x))./(1+((x-position)./(0.5.*width)).^2);