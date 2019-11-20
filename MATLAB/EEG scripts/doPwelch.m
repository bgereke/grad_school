function [P,F] = doPwelch(D,Fs,nfft)
%
% Calculate powerspectra using Welch method (50% overlap windows; Hanning
% taper)
%
% nfft determines spectral smooting: spectral smooting decreases with
% increasing nfft. Use n^2 values (512, 1024, etc)
% Default nfft = 4096
%
%
% IMPORTANT: Remember to set sampling frequency Fs 
%
% Ole Jensen, May 10, 2006

%Fs = 2034;
%Fs = 1893;

if nargin <= 2
    nfft = 4096;
end
[P,F] = pwelch(D,nfft,nfft/2,nfft,Fs);

%
% Plotting examples:
%
% figure(1)
% plot(F,P); xlabel('frequency (Hz)'); ylabel('power');
% xlim([0 100])
% figure(2)
% plot(F,log(P)); xlabel('frequency (Hz)'); ylabel('log(power)');
% xlim([0 100])
