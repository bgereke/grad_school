function [Pxx,f] = psdd(x,Fs);

Pxx  =  abs(fft(x));
f = (0:length(x)-1)*Fs/length(x);
f = f(1:floor(length(f)/2));
Pxx = Pxx(1:floor(length(Pxx)/2));
if (nargout == 0) 
   plot(f,10*log10(Pxx)), grid on
   xlabel('Frequency'), ylabel('Power Spectrum Magnitude (dB)');
end 
