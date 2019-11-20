function [eeg_epochs_grouped, sPhase, sPhase2] = spike_triggered_windows_NLX_gamma40_phasevec(eeg_ind_spikes, EEG, f1, f2, Fs)

%  eeg_ind_spikes are the spike time stamps, converted to align with the EEG points
%  EEG is a single EEG trace
%  Fs is sampling frequency

%  eeg_epochs_grouped contains the spike-triggered windows for STAs
%  sPhase contains the phase estimates, one for each spike time stamp

    eeg_epochs_grouped = [];
    
    phase = [];
    
    sPhase = [];
    
    frequency_for_phase_vec = (f1+f2)/2;
    
    %B1 = phasevec_t(frequency_for_phase_vec,detrend(EEG),Fs,9);
    %B1 = phasevec_t(frequency_for_phase_vec,detrend(EEG),Fs,5);
    %B1 = phasevec_t(frequency_for_phase_vec,detrend(EEG),Fs,7);
    
    %phase = angle(B1);
    
    
    
    phase2 = [];
    
    sPhase2 = [];
    
    EEG_gamma_filter = fftbandpass(EEG, Fs, f1 -2, f1, f2, f2+2);
    %EEG_gamma_filter = bandpass(EEG, f1, f2, Fs);
    DTAS = hilbert(EEG_gamma_filter);
    phase2 = angle(DTAS);
    
    
    for k = 1:length(eeg_ind_spikes)

        
        eeg_epoch_to_insert = EEG(eeg_ind_spikes(k) - 400: eeg_ind_spikes(k) + 400);
        eeg_epochs_grouped = [eeg_epochs_grouped, eeg_epoch_to_insert];
    
        
        %sPhase = [sPhase,phase(eeg_ind_spikes(k))];
        
        sPhase2 = [sPhase2,phase2(eeg_ind_spikes(k))];
        
    end