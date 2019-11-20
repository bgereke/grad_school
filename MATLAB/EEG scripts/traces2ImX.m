function [ImX] = traces2ImX(x,y,freqVec,Fs,width)
% function [WPLI,timeVec,freqVec] = traces2ImC(S,freqVec,Fs,width);
%
% Calculates the time frequency representation of the imaginary component
% of the cross-spectrum between two signals using a Morlet wavelet method. 
% see Lachaux et al. Neurophysiol Clin 2002 ; 32: 157-74
%
% Input
% -----
% x,y    : signals (time x 1)  
% freqVec    : frequencies over which to calculate TF energy        
% Fs   : sampling frequency
% cm: number of cycles in morlet wavelet (> 5 advisable)  
% ci: number of cycles in cross-spectral integration window (>= cm advisable, <= 6 for 100ms resolution)
%
% Output
% ------
% ImX    : Imaginary component of cross-spectrum(frequency x time)
%     
%------------------------------------------------------------------------
% Ole Jensen, Brain Resarch Unit, Low Temperature Laboratory,
% Helsinki University of Technology, 02015 HUT, Finland,
% Report bugs to ojensen@neuro.hut.fi
%------------------------------------------------------------------------

%    Copyright (C) 2000 by Ole Jensen 
%    


x = x'; y = y';
Wx = zeros(length(freqVec),size(x,2)); 
Wy = zeros(length(freqVec),size(x,2));
ImX = zeros(length(freqVec),size(x,2));

for j=1:length(freqVec)
    Wx(j,:) = energyvec(freqVec(j),detrend(x),Fs,width);
    Wy(j,:) = energyvec(freqVec(j),detrend(y),Fs,width);
end

% for j=1:length(freqVec)
%     ImX(j,:) = xspec(freqVec(j),Wx(j,:),Wy(j,:),Fs,width);
% end

ImX = imag(Wx.*conj(Wy));
% ImX = imag(ImX);


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
y = convfft(s,m);
%y = conv(s,m);
%y = (2*abs(y)/Fs).^2;
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
A = (st*sqrt(pi))^(-0.5);

y = A*exp(-t.^2/(2*st^2)).*exp(i*2*pi*f.*t);


function X = xspec(f,Wx,Wy,Fs,width)
% function X = xspec(f,Wx,Wy,Fs,width)
%
% Return a vector containing the cross-spectrum as a
% function of time for frequency f. The cross-spectrum
% is calculated by integrating over a window around time t
% set by width. 
% Wx,Wy : Morlet-transformed signals
% Fs: sampling frequency
% width : width of integtration window (>= 5 suggested).
%
% Brian Gereke, February 2015

dt = 1/Fs;
sf = f/width;
st = 1/(2*pi*sf);

t=-3.5*st:dt:3.5*st;
g = gaussian(f,t,width);
X = Wx.*conj(Wy);
% X = convfft(X,g);
% %y = conv(s,m);
% %y = (2*abs(y)/Fs).^2;
% X = X(ceil(length(g)/2):length(X)-floor(length(g)/2));

function g = gaussian(f,t,width)
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
A = (st*sqrt(pi))^(-0.5);

g = A*exp(-t.^2/(2*st^2));

