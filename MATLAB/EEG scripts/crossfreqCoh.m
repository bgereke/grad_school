function [CRF,freqVec1,freqVec2] = crossfreqCoh(S1,S2,freqVec,Fs);
% function [CRF,f] = crossfreqCoh(S1,S2,freqVec,Fs);
%
% [CRF,freqVec1,freqVec2] = crossfreqCoh2(Data,Data,2:2:200,2000);
% imagesc(freqVec1,freqVec2,CRF);axis xy;
% xlim([0 100])


width = 5;

%nfft = 2048*2; % Choose about 4*Fs
nfft = 2048*4; % Choose about 4*Fs
% nfft = 2048/2; % Choose about 4*Fs


freqVec2 = freqVec';
for k=1:length(freqVec)
    %fprintf('%d ',freqVec(k))
    E1 = energyvec(freqVec(k),S1',Fs,width);
    [CRF(k,:),freqVec1] = mscohere(E1',S2,hanning(nfft),nfft/2,nfft,Fs);
end
fprintf('\n');

CRF = CRF(:,1:floor(end/2));
freqVec1 = freqVec1(1:floor(end/2));

function y = energyvec(f,s,Fs,width)
% function y = energyvec(f,s,Fs,width)
%
% Return a vector containing the energy as a
% function of time for frequency f. The energy
% is calculated using Morlet's wavelets. 
% s : signal
% Fs: sampling frequency
% width : width of Morlet wavelet (>= 5 suggested).
%
% 

dt = 1/Fs;
sf = f/width;
st = 1/(2*pi*sf);

t=-3.5*st:dt:3.5*st;
m = morlet(f,t,width);
y = conv(s,m');
y = (2*abs(y)/Fs).^2;
y = y(ceil(length(m)/2):length(y)-floor(length(m)/2));

function y = morlet(f,t,width)
% function y = morlet(f,t,width)
% 
% Morlet's wavelet for frequency f and time t. 
% The wavelet will be normalized so the total energy is 1.
% width defines the ``width'' of the wavelet. 
% A value >= 5 is suggested.
%
% Ref: Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997)
%
% See also: PHASEGRAM, PHASEVEC, WAVEGRAM, ENERGY 
%
% Ole Jensen, August 1998 


sf = f/width;
st = 1/(2*pi*sf);
A = 1/(st*sqrt(2*pi));

y = A*exp(-t.^2/(2*st^2)).*exp(i*2*pi*f.*t);

