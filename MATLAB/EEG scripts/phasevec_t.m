
function y = phasevec_t(f,s,Fs,width)
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
%for k=1:size(s,2) 
    y = conv(s,m');
%end
l = find(abs(y) == 0);
y(l) = 1;
y = y./abs(y);
y(l) = 0;
y = y(ceil(length(m)/2):length(y)-floor(length(m)/2),:);


