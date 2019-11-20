function [eeg_epochs_grouped, eeg_epoch_index2] = find_spike_times_rs16_STA(TS_spikes, ts_eeg, EEG, Fs)

% TS_spikes is the time stamp for the spikes, output of the 'loadSpikes'
% function.
% ts_eeg is the time stamp for the EEG recordings (CSC file), one time
% stamp for every 512 points.  This variable is part of the output of the
% 'load_eeg' function.
% EEG is the CSC file containing the EEG recording to which the spikes will be aligned 

eeg_epochs_grouped = [];

eeg_epoch_index2 = [];

ts2 = ts_eeg / 1000000;  % converts spike time stamp from microseconds to seconds

%eliminate spikes that occur earlier than the first time stamp of the EEG
% TS_spikes_index = find (TS_spikes > ts2(1));  
% TS2_spikes = TS_spikes (TS_spikes_index(1) : length(TS_spikes) );
TS2_spikes = TS_spikes(TS_spikes >= ts2(1));

%EEG_theta_filter = bandpass(EEG, 6, 10, Fs );
%DTAS = hilbert(EEG_theta_filter);
%phase = angle(DTAS);

nSpikes = length(TS2_spikes);
%sPhase = zeros(nSpikes,1);

for k = 1:length(TS2_spikes)
    n = TS2_spikes(k);
%     [eeg_index_boundary] = find(n < ts2 & ts2 > n + 0.2705);

    % Find closest eeg timestamp to the current spike timestamp
    tDiff = (ts2-n).^2;
    [minDist,eeg_index_boundary] = min(tDiff);
    if length(eeg_index_boundary) > 1
        eeg_index_boundary = eeg_index_boundary(1);
    end
    
    % Time difference between eeg timestamp and spike timestamp
    diff = ts2(eeg_index_boundary) - n;
    
    if diff == 0 % Hit timestamp right on
        eeg_epoch_index = 512*(eeg_index_boundary-1)+1;
    elseif diff < 0 % Spike timestamp is in next 512 point data segment
        offset = round(abs(diff)*Fs);
        eeg_epoch_index = 512*(eeg_index_boundary-1) + offset;
    else % Spike timestamp is in last 512 point data segment
        offset = round(diff*Fs);
        eeg_epoch_index = 512*(eeg_index_boundary-1) - offset;
    end
    
%     samples_conv = diff * Fs;
%     eeg_epoch_index = (512 * eeg_index_boundary) + samples_conv;
    
    if eeg_epoch_index < (Fs*.2)
        
        eeg_epochs_grouped = eeg_epochs_grouped;
        
        eeg_epoch_index2 = eeg_epoch_index2;
        
    elseif eeg_epoch_index > (length(EEG)-(Fs*.2))
        
        eeg_epochs_grouped = eeg_epochs_grouped;
        
        eeg_epoch_index2 = eeg_epoch_index2;
    
    else 
        
        eeg_epoch_to_insert = EEG(eeg_epoch_index - (Fs*.2): eeg_epoch_index + (Fs*.2));
        eeg_epochs_grouped = [eeg_epochs_grouped, eeg_epoch_to_insert];
        
        
        eeg_epoch_index2 = [eeg_epoch_index2, eeg_epoch_index];
                
    end
    
    %sPhase(k) = phase(eeg_epoch_index);
    
    
    
end


    
   

