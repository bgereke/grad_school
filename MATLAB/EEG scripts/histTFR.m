function [TFRnorm, TFRcounts] = histTFR(TFR,freq,powbins)

TFRnorm = TFR./repmat(max(TFR,[],2),1,size(TFR,2));
TFRcounts = histc(TFRnorm,powbins,2);
% figure
% imagesc(powbins,freq,TFRcounts,[0 0.005*max(max(TFRcounts))]); axis xy
% xlabel('normalized power');ylabel('frequency (Hz)');title('Power-Frequency Histogram')
% colorbar
% figure
% plot(freq,TFRcounts(:,1)/2000);xlabel('frequency (Hz)');ylabel('time (sec)'); % local minima denote frequencies of interest
% figure
% plot(freq,max(TFR,[],2));xlabel('frequency (Hz)');ylabel('max power');
% v = TFRcounts(max(freq),1:length(powbins));
% figure
% contour(powbins,freq,TFRcounts,v);