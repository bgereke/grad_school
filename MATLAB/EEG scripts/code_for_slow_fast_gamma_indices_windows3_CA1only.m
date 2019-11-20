function [peak_ind_unique,non_overlap_peak_ind_unique,start2,stop2,non_overlap_gamma_windows_CA1] = code_for_slow_fast_gamma_indices_windows3_CA1only(EEG,EEG2, f1, f2, Fs);

%EEG = ch3_rev_begin1;

%f1 = 100;    %  frequency band of interest
%f2 = 140;

[TFR] = TFR_frequency_band(EEG,Fs,5, f1, f2);
[wpli] = traces2WPLI(EEG2,EEG,f1:2:f2,Fs,6);
%TFR = TFR./repmat(max(TFR,[],2),1,size(TFR,2));
TFR_z = zscore(TFR);
high = find(TFR_z > 1 & mean(wpli,1)>0.8);
bp = fftbandpass(EEG,Fs,f1-2,f1,f2,f2+2);
%bp = bandpass(EEG,f1,f2,Fs);
detect_length = round(0.16 * Fs); %   length of the window of interest (time * sampling frequency)
detect_start = round(0.08 * Fs);  %   start of the window of interest
kc = 0;
for k = 1:size(high,2)              %   for all samples with a high power
start = high(k) - detect_start;    %   determine start of window of interest
stop = start + detect_length;               %   determine stop of window of interest
if 0 < start & length(bp) >= stop        %   window has to be inside the whole EEG
kc = kc + 1;                            %   counter for windows of interest inside the whole EEG
[max_val max_ind]= max(bp(start:stop));  %   determine the index of the maximal value in the bandpassed EEG
peak_ind(kc) = start + max_ind - 1;      %   add the offset, in peak_ind are now all the sample indices with high power
end
end
p = 0;
counter = 0;
for k = 1:size(peak_ind,2)
if p < peak_ind(k)                      %   we don't want duplicates
p = peak_ind(k);
start = p - round(0.2 * Fs);      %   start of the windows that become averaged
stop = start + round(0.4 * Fs);   %   stop of the windows that become averaged
%start = p - (round(Fs/10));      %   start of the windows that become averaged
%stop = p + (round(Fs/10));   %   stop of the windows that become averaged
if 0 < start & length(EEG) >= stop %   windows have to be inside the whole EEG
counter = counter + 1;
peak_ind_unique(counter) = p;   %   here are the indices with eliminiated duplicates
%gamma_windows_25_140_z2(counter,:) = EEG(start:stop);
%slow_gamma_windows_ca1(counter,:) = egf_27110711(start:stop);
%fast_gamma_windows_ec(counter,:) = t_ec.EEG_rev_uvolts(start:stop);
%fast_gamma_windows_ca3(counter,:) = t_ca3.EEG_rev_uvolts(start:stop);
end
end
end


non_overlap_peak_ind_unique = [];
for n=2:length(peak_ind_unique)
if peak_ind_unique(n-1) + (round(Fs/10)) < peak_ind_unique(n);
non_overlap_peak_ind_unique = [non_overlap_peak_ind_unique,peak_ind_unique(n-1)];
else
non_overlap_peak_ind_unique = non_overlap_peak_ind_unique;
end
end

p = 0;
counter = 0;
for k = 1:size(non_overlap_peak_ind_unique,2)
if p < non_overlap_peak_ind_unique(k)                      %   we don't want duplicates
p = non_overlap_peak_ind_unique(k);
start = p - round(0.2 * Fs);
stop = start + round(0.4 * Fs);
if 0 < start & length(EEG) >= stop
counter = counter + 1;
non_overlap_gamma_windows_CA1(counter,:) = EEG(start:stop);
%non_overlap_gamma_windows_EC(counter,:) = EEG_EC(start:stop);
%non_overlap_gamma_windows_CA3(counter,:) = EEG_CA3(start:stop);

start2(counter) = start;
stop2(counter) = stop;

end
end
end


% plot((1:size(non_overlap_gamma_windows_CA1,2))/Fs*1000,mean(non_overlap_gamma_windows_CA1));

%hold on;plot(mean(non_overlap_gamma_windows_CA3),'g');




