function S = LoadSpikes_klusta(session,numShanks)

S = cell(1);
cellnum = 1;

for shank = 1:numShanks
    %load spike times and cluster information from .kwik files
    filename = strcat(session(1:end-7), 'S', num2str(shank), 'CSCs.kwik');
    spkidx=hdf5read(filename, '/channel_groups/0/spikes/time_samples');
    clnum=hdf5read(filename, '/channel_groups/0/spikes/clusters/main');
    cluList=unique(clnum);
    %iterate through clusters
    for k = 1:length(cluList)
        %only take 'good' clusters
        clGroup = hdf5read(filename,['/channel_groups/0/clusters/main/' num2str(cluList(k))],'cluster_group');
        if clGroup == 2
            %store cluster spike indices
            S{cellnum} = spkidx(clnum==cluList(k));
            cellnum = cellnum+1;
        end
    end
end