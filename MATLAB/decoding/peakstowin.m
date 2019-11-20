function [times] = peakstowin(peakind,eeg,eegTS,winsize)

%use for getting the windows around the peaks, when the eeg has
%been downsampled (by eegendoftrackcutter etc) and the windows may be 
%straddling cutoff times

times(:,1) = peakind - round(winsize*2000/2);
times(:,2) = peakind + round(winsize*2000/2);

%make sure all the points are within the boundaries
for i = 1:size(times,1)
    if times(i,1) < 1;
        times(i,1) = 1;
    end
    if times(i,2) > length(eeg)
        times(i,2) = length(eeg);
    end
end

%find straddled windows and truncate them
jumps = find(diff(eegTS)>(1/500));
for i = 1:size(times,1)
    ind = find(jumps>=times(i,1) & jumps <=times(i,2));
    if length(ind)==1
        if jumps(ind) - times(i,1) > (winsize*2000/2)
            times(i,2) = jumps(ind) - 1;
        else
            times(i,1) = jumps(ind) + 1;
        end
    elseif length(ind)>1
        times(i,:) = nan;
    end
end
times(isnan(times(:,1)),:) = [];
times = eegTS(times);