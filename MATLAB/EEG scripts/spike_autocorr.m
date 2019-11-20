[TS] = loadSpikes('C:\Data\mouse11\2013-10-17_14-30-03\begin1\TT1_3.t');

Fs = 2000;
intervals = diff(TS);
no_spikes = floor(intervals*Fs);
train = [1];

for i=1:length(no_spikes)
    train = [train zeros(1,no_spikes(i)) 1];
end

bin = 0.001; % in sec
indexed_bin = bin*Fs;

for ii=1:floor(length(train)/indexed_bin)
        binned_train(ii) = sum(train((indexed_bin*(ii-1)+1):indexed_bin*ii));
end

[xc, lags] = xcorr(binned_train,floor(0.3/bin),'coeff');
plot(lags*bin, xc)

xlabel('time lag (sec)')
ylabel('normalized autocorrelation')
title('mouse11 101713 begin1 TT1_3 autocorrelation')
xlim([-0.3 0.3])
ylim([0 0.07])
