% Finds the position of the spikes
function [pos] = spikePos(ts,phase,post)
N = length(ts);
pos = zeros(N,1);

for ii = 1:N
    tdiff = (post-ts(ii)).^2;
    [~,ind] = min(tdiff);
    pos(ii) = phase(ind(1));
end