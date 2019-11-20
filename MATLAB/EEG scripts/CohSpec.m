function [Coh,freqVec] = CohSpec(S1,S2,freqVec,Fs);

width = 5;

%nfft = 2048*2; % Choose about 4*Fs
%nfft = 2048*4; % Choose about 4*Fs
nfft = 2048; % Choose about 4*Fs


[Coh,freqVec] = mscohere(S1,S2,hanning(nfft),nfft/2,nfft,Fs);


