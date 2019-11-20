function [B,timeVec,freqVec] = traces2PLV_1trial(S1,S2,freqVec,Fs,width);
% function [B,FREQVEC,TIMEVEC] = TRACES2PLV(S1,S2,FREQVEC,FS,WIDTH);
% 
% Calculats the 'phase-locking statistics' over multiple sets of trials using 
% Morlet wavelets.
%
% S1       : time x trials  (signal 1)    
% S2       : time x trials  (signal 2)
% FREQVEC  : vector of frequency points   
% TIMEVEC  : vector  time points
% FS       : sampling frequency
% WIDTH    : width of Morlet wavelet (>= 5 suggested).
% B        : Phase synchrony (0 to 1), frequnecy  x time
% BSTAT    : Phase-locking statistics, frequency x time
% BPLF     : Stimulus locked phase-locking
%
% Lachaux JP, Rodriguez E, Martinerie J, Varela FJ.
% Measuring phase synchrony in brain signals.
% Hum Brain Mapp. 1999;8(4):194-208.
%
%------------------------------------------------------------------------
% Ole Jensen, Brain Resarch Unit, Low Temperature Laboratory,
% Helsinki University of Technology, 02015 HUT, Finland,
% Report bugs to ojensen@neuro.hut.fi
%------------------------------------------------------------------------

%    Copyright (C) 2000 by Ole Jensen 
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You can find a copy of the GNU General Public License
%    along with this package (4DToolbox); if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

timeVec = (1:length(S1))/Fs;  
B = zeros(length(freqVec),length(S1));

%fprintf('Frequency (Hz):\n');
for j=1:length(freqVec) 
%    fprintf('%d ',freqVec(j));
   B1 = phasevec(freqVec(j),detrend(S1),Fs,width);
   B2 = phasevec(freqVec(j),detrend(S2),Fs,width);

      B(j,:) = (B1./B2)';                         

end


function y = phasevec(f,s,Fs,width)
% function y = phasevec(f,s,Fs,width)
%
% Return a the phase as a function of time for frequency f.
% The phase is calculated using Morlet's wavelets.
%
% Fs: sampling frequency
% width : width of Morlet wavelet (>= 5 suggested).
%
% Ref: Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997)

dt = 1/Fs;
sf = f/width;
st = 1/(2*pi*sf);
t=-3.5*st:dt:3.5*st;
m = morlet(f,t,width);
for k=1:size(s,2) 
    y(:,k) = convfft(s(:,k),m);
end
l = find(abs(y) == 0);
y(l) = 1;
y = y./abs(y);
y(l) = 0;
y = y(ceil(length(m)/2):size(y,1)-floor(length(m)/2),:);

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
%
% Ole Jensen, August 1998

sf = f/width;
st = 1/(2*pi*sf);
A = 1/sqrt(st*sqrt(pi));
y = A*exp(-t.^2/(2*st^2)).*exp(i*2*pi*f.*t);

