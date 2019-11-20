function [ifdgram,sonogram]=ifdv(s,SAMPLING,N,OVERLAP,sigma,ZOOMT,ZOOMF,TL,FL)

%Remapped sonograms, as described in
%Gardner & Magnasco PNAS 2006
%
%For lack of a better name, at the moment, "ifdgram" is the new sonogram
% 
% [s] Signal to be analyzed.
% [SAMPLING] Sampling rate of the signal. 
% [N] Number of filters in the filter-bank.
%------
% [OVERLAP] Number of samples to overlap in each successive window;
%
% The number of time points in the final image is proportional to
% length(s)/(N-OVERLAP);
% Higher overlap will result in sharper lines.
%------
% [sigma] is the temporal resolution of the analysis in milliseconds.
% N should be larger than 5*SAMPLING*(sigma/1000);
%
% Choose sigma small to represent sound in a time-like fashion - as a
% series of clicks, or sigma large, to represent sound in a frequency-like
% fasion, as a series of tones. For most signals of interest, intermediate
% values of sigma are best.
%-------
% [ZOOMT] sets the temporal resolution of the final image typically ZOOMT=1.
%-------
%[ZOOMF] sets the resultion of final image in Frequency. Res=ZOOMF*(N/2)
% In contrast to ZOOMT, it is typically useful to set ZOOMF great than one,
%------
%
%[TL] temporal locking window in pixels
%[FL] frequency locking window in pixels
%
%When the remapping moves a pixel by more than TL, or FL, that pixel
%acquires zero weight. For discussion of the locking window, see
%Gardner & Magnasco, J. Acoust. Soc. Am. 2005
%
%When these parameters are small (order 1) "stray" points are removed and the lines are sharpened.
%If these parameters are too small, lines become too thin, and appear
%discontinuous.
%
%Typical parameters: SAMPLING 44100, N=1024, OVERLAP=1010, 
%sigma=2, ZOOMT=1, ZOOMF=3; FL=5;TL=5
%---
% Implementation note:
%The best results will come from calculating an ifdgram for many values of
%sigma, (.5:.1:3.5) for example, then combining by multiplying together images
%with neighboring values of sigma, and adding them all together.
%Rational for this is given in Gardner & Magnasco PNAS 2006.
%---------------
%example: Compute an ifdgram of 100ms of white noise
%s=rand(4000,1)-0.5;
%[ifdgram,sonogram]=ifdv(s,44100,1024,1020,1,1,1,2,2);
%colormap(hot)
%imagesc(log(ifdgram+3));
%-------
%Comment: As in the previous example, log scaling of the ifdgram,
%may be optimal for most sounds.

FACTOR=double(N-OVERLAP);
t=-N/2+1:N/2;
%Gaussian and first derivative as windows.

sigma=(sigma/1000)*SAMPLING;
w=exp(-(t/sigma).^2);
dw=(w).*((t)/(sigma^2))*-2;

q=specgram(s,N,[],w,OVERLAP)+eps; %gaussian windowed spectrogram
q2=specgram(s,N,[],dw,OVERLAP)+eps; %deriv gaussian windowed spectrogram
[F,T]=size(q);
dx=(q2./q)/(2*pi); %displacement according to the remapping algorithm
sonogram=q;
fo=zeros(size(q));
to=zeros(size(q));
for(k=1:T) to(:,k)=k; end
for(n=1:F) fo(n,:)=(.5*n-1)/(F-1); end

f_est= ((fo-imag(dx))*N)+1; %calculate frequency displacement factors
t_est= to-(pi*sigma*sigma)*real(dx)/FACTOR; %calculate temporal displacement

tref=ZOOMT*repmat(1:T,F,1);
fref=ZOOMF*repmat(1:F,T,1)';
F=F*ZOOMF;T=T*ZOOMT; %rescale the image by the resolution stretching factors.
f_e=round(ZOOMF*f_est);
t_e=round(ZOOMT*t_est);

q(f_e<1 | f_e>F)=0; %set to zero the points that are mapped out of the image
q(t_e<1 | t_e>T)=0;

f_e(f_e<1)=1; f_e(f_e>F)=F; t_e(t_e<1)=1; t_e(t_e>T)=T;%keep all points in the image
q(abs(f_e-fref)>FL)=0; q(abs(t_e-tref)>TL)=0; %remove "stray points"
[a,b]=size(f_e);n=a*b;
newq=accumarray([f_e(1:n)' t_e(1:n)'],abs(q(1:n)')); %build the remapped sonogram

%the Frequency axis of the final image is zero to Nyquist freqeuncy, that is
%(SAMPLING/2) Hz.

ifdgram=abs(flipdim(newq,1)); %flip image, so low frequency is at the bottom.
sonogram=abs(flipdim(sonogram,1)); %this is the standard sonogram



