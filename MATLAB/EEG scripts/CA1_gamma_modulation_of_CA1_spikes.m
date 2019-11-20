function [eeg_epoch_index2, eeg_epoch_index2_gamma_epochs_unique, sPhase2_CA1_cell_ca1_gamma_epochs, phaseBin_sPhase2_CA1_cell_ca1_gamma_epochs, phaseBin2_sPhase2_CA1_cell_ca1_gamma_epochs] = CA1_gamma_modulation_of_CA1_spikes2(TS_spikes, ts_eeg, eeg_CA1, peak_ind_unique, f1, f2, Fs);


[eeg_epoch_index2] = find_spike_times_rs13_no_phase(TS_spikes, ts_eeg, eeg_CA1, Fs);
eeg_epoch_index2_gamma_epochs = [];

p = 0;
counter = 0;
counter2 = 0;
for k = 1:length(peak_ind_unique)

if p < peak_ind_unique(k)                      %   we don't want duplicates
p = peak_ind_unique(k);
start = p - round(0.2 * Fs);      %   start of the windows that become averaged
stop = start + round(0.4 * Fs);   %   stop of the windows that become averaged
if 0 < start && length(eeg_CA1) >= stop %   windows have to be inside the whole EEG
counter = counter + 1;
%peak_ind_unique(counter) = p;   %   here are the indices with eliminiated duplicates
%gamma_windows_hipp(counter,:) = ch3_rev_begin1(start:stop);
%gamma_windows_EC(counter,:) = ch12_rev_begin1(start:stop);
for m = 1:length(eeg_epoch_index2)
if eeg_epoch_index2(m) > start && eeg_epoch_index2(m) < stop 
counter2 = counter2 + 1;
eeg_epoch_index2_gamma_epochs(counter2) = eeg_epoch_index2(m);
end
end
end
end
end

if ~isempty(eeg_epoch_index2_gamma_epochs)
    eeg_epoch_index2_gamma_epochs_unique = unique(eeg_epoch_index2_gamma_epochs);
    
    [eeg_epochs_grouped_ca1_gamma_epochs, sPhase_CA1_cell_ca1_gamma_epochs, sPhase2_CA1_cell_ca1_gamma_epochs] = spike_triggered_windows_NLX_gamma40_phasevec(eeg_epoch_index2_gamma_epochs_unique, eeg_CA1, f1, f2, Fs);
    
    %[phaseBin_sPhase_CA1_cell_ca1_gamma_epochs, phaseBin2_sPhase_CA1_cell_ca1_gamma_epochs] = make_spike_time_histogram_rep(sPhase_CA1_cell_ca1_gamma_epochs);
    [phaseBin_sPhase2_CA1_cell_ca1_gamma_epochs, phaseBin2_sPhase2_CA1_cell_ca1_gamma_epochs] = make_spike_time_histogram_rep(sPhase2_CA1_cell_ca1_gamma_epochs);
else
    eeg_epoch_index2, eeg_epoch_index2_gamma_epochs_unique = [];
    sPhase2_CA1_cell_ca1_gamma_epochs = [];
    phaseBin_sPhase2_CA1_cell_ca1_gamma_epochs = [];
    phaseBin2_sPhase2_CA1_cell_ca1_gamma_epochs = [];
end
%done