% sPhase is part of the output of the find_spike_times_rs3 function and
% contains the extracted phases at the particular times at which spikes occurred



function [phaseBin, phaseBin2] = make_spike_time_histogram_rep(sPhase)

% Convert sPhase to degrees
sPhase2 = (sPhase+pi) * 360/(2*pi);
sPhase3 = ((sPhase+pi) * 360/(2*pi)) + 360;

start = 0;
stop = 30;
phaseBin = zeros(12,1);
for ii=1:12
ind = find(sPhase2 >= start & sPhase2 < stop);
phaseBin(ii) = length(ind);
start = start + 30;
stop = stop + 30;
end


start2 = 360;
stop2 = 390;
phaseBin2 = zeros(12,1);
for ii = 1:12
ind2 = find(sPhase3 >= start2 & sPhase3 < stop2);
phaseBin2(ii) = length(ind2);
start2 = start2 + 30;
stop2 = stop2 + 30;
end