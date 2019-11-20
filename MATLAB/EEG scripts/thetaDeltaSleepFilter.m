function [pass_times, pass_idx] = thetaDeltaFilter(eeg, tt, Fs,tsp, td_threshold)

%make the filtered eeg traces
winsize=2*Fs;  %looks at 2 second, non-overlapping windows
theta_eeg = fftbandpass(eeg,Fs, 4,5,10,11);
delta_eeg = fftbandpass(eeg,Fs,1,2,4,5);
td_total = sum(periodogram(theta_eeg,[],[],Fs))/sum(periodogram(delta_eeg,[],[],Fs));

pass_times=[];
pass_idx = [];
td_max = nan;
td_min = nan;
td_mean=0;
trialend=min(length(eeg),length(tt))-winsize;

for aa=1:winsize:trialend
	
	win_theta_eeg=theta_eeg(aa:aa+winsize-1);
	win_delta_eeg = delta_eeg(aa:aa+winsize-1);
	td_rat = sum(periodogram(win_theta_eeg,[],[],Fs))/sum(periodogram(win_delta_eeg,[],[],Fs));
	td_max=nanmax(td_rat,td_max);
	td_min=nanmin(td_rat,td_min);
	td_mean=td_mean + td_rat/(length(eeg)/winsize);
	
	win_amp = mean(abs(eeg(aa:aa+winsize-1)));
    
	if td_rat > td_threshold 
        pass_times = [pass_times;tt(aa), tt(aa+winsize-1)];
        pass_idx = [pass_idx;aa, aa+winsize-1];        
    end
    
    if pass_idx(end,1) == pass_idx(end-1,2)+1
       pass_idx(end-1,2) = pass_idx(end,2);
       pass_times(end-1,2) = pass_idx(end,2);
       pass_idx(end,:) = [];
       pass_times(end,:) = [];
    end
end
% if(isempty(pass_times))
% 	tsp_filt=[];
% else
% 	[tsp_filt, ~]= get_spike_times_within_specific_eeg_windows(pass_times(:,1),pass_times(:,2),tsp);
% end