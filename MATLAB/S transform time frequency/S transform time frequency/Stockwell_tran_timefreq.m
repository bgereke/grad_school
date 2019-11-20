function [  ] = Stockwell_tran_timefreq(a,minfreq,maxfreq,samplingrate)
%Time frequency distribution using S-transform (Stockwell transform)
%'a' is the input signal or time series
% 'minfreq' is the minimum frequency of the signal
% 'maxfreq' is the maximum frequency of the signal (sampling rate/2)
% 'samplingrate' is the sampling rate of the signal
% The Y axis is frequency going from minfrequency to maxfrequency
% The X axis is time/ Sample
%For any further queries contact aditsundar@gmail.com

[st1,t,f] = st(a,minfreq,maxfreq,samplingrate,1) ;

figure(1)
plot(a)
z=abs(st1);
z= imcomplement(z);
imtool(z)    


end

