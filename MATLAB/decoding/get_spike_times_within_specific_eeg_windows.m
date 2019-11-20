function  [episode_spike_time, episode_spike_time_index, spksfgz] = get_spike_times_within_specific_eeg_windows(startTime, stopTime, spike_time, wt, sgz, fgz)

sampCounter = 0;

for jj = 1:length(startTime)
        temp = spike_time;
        ind = find(temp >= startTime(jj) & temp <= stopTime(jj));
        N = length(ind);
        episode_spike_time(sampCounter+1:sampCounter+N) = wt(jj)*ones(N,1);
        episode_spike_time_index(sampCounter+1:sampCounter+N) = ind;
        spksfgz(sampCounter+1:sampCounter+N,1) = sgz(jj)*ones(1,N);
        spksfgz(sampCounter+1:sampCounter+N,2) = fgz(jj)*ones(1,N);
        sampCounter = sampCounter + N;
end
    