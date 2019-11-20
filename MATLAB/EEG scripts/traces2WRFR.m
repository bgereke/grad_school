function [wrfr] = traces2WRFR(x,y,freqVec,Fs,width);
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
wrfr = zeros(length(freqVec),size(x,2));

for j=1:length(freqVec)
    Wx = energyvec(freqVec(j),detrend(x),Fs,width);
    Wy = energyvec(freqVec(j),detrend(y),Fs,width);
    X = Wx.*conj(Wy);
    R = abs(X);
    ImX = imag(X);
    sgn = sign(mean(ImX));
    ImX(sign(ImX)~=sgn) = 0;
    ImX = abs(ImX);
    wrfr(j,:) = waverfr(freqVec(j),R,ImX,Fs,50);
end

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

y = A*exp(-t.^2/(2*st^2)).*exp(1i*2*pi*f.*t);


function wrfr = waverfr(f,ImX,Fs,width)
% function wpli = wpli(f,Wx,Wy,Fs,width)
%
% Return a vector containing the wpli as a
% function of time for frequency f. The wpli
% is calculated by integrating over a gaussian window around time t
% set by width. 
% ImX: imaginary part of the cross spectrum (frequency x time)
% Fs: sampling frequency
% width : width of integtration window (>= 5 suggested).
%
% Brian Gereke, February 2015

dt = 1/Fs;
sf = f/width;
st = 1/(2*pi*sf);

t=-3.5*st:dt:3.5*st;
g = gaussian(f,t,width);
num = convfft(R.*ImX,g);
den = convfft(ImX,g);
% ossq = convfft(ImX.^2,g);
num = num(ceil(length(g)/2):length(num)-floor(length(g)/2));
den = den(ceil(length(g)/2):length(den)-floor(length(g)/2));
% ossq = ossq(ceil(length(g)/2):length(ossq)-floor(length(g)/2));
wrfr = num./den;

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

