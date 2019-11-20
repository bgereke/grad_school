file = 'C:\Data\mouse11\2013-10-17_14-30-03\begin1\CSC1.ncs';
[TS] = loadSpikes('C:\Data\mouse11\2013-10-17_14-30-03\begin1\TT1_4.t')
load_eeg5;
ch_X = samples;  
%load_eeg_4_retrieve_header;
header
%you'll see the ADBitVolts conversion factor.... e.g., ADBitVolts 0.00644
ch_X = ADBitVolts .* ch_X;  %convert to volts, get ADBitVolts from load_eeg_4_retrieve_header
ch_X = ch_X .* 1000000;

f1 =65; f2 = 100;

[peak_ind_unique,non_overlap_peak_ind_unique,start2,stop2,non_overlap_gamma_windows_CA1] = code_for_slow_fast_gamma_indices_windows3_CA1only(ch_X, f1, f2, 2000);
[eeg_epoch_index2, eeg_epoch_index2_gamma_epochs_unique, sPhase2_CA1_cell_ca1_gamma_epochs, phaseBin_sPhase2_CA1_cell_ca1_gamma_epochs, phaseBin2_sPhase2_CA1_cell_ca1_gamma_epochs] = CA1_gamma_modulation_of_CA1_spikes(TS, ts, ch_X, peak_ind_unique, f1, f2, 2000);

% [eeg_epoch_index2] = find_spike_times_rs13_no_phase(TS, ts, ch_X, 2000);
% [eeg_epochs_grouped_ca1_gamma_epochs, sPhase_CA1_cell_ca1_gamma_epochs, sPhase2_CA1_cell_ca1_gamma_epochs] = spike_triggered_windows_NLX_gamma40_phasevec(eeg_epoch_index2, ch_X, f1, f2, 2000)
% [phaseBin_sPhase2_CA1_cell_ca1_gamma_epochs, phaseBin2_sPhase2_CA1_cell_ca1_gamma_epochs] = make_spike_time_histogram_rep(sPhase2_CA1_cell_ca1_gamma_epochs);

bar([0:30:690],[phaseBin_sPhase2_CA1_cell_ca1_gamma_epochs;phaseBin2_sPhase2_CA1_cell_ca1_gamma_epochs],1)
xlabel('phase (deg)')
ylabel('spike count')
title('mouse11 082513 begin1 CSC1 TT1_4 fast gamma histogram')


