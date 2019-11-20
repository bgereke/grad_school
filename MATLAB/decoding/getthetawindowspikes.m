function [result] = getthetawindowspikes(fields,thetawindow,minspikes)

%create matrix of spike time and the center of mass (com) nesxt to it then
%sort by spike time


numwindows = size(thetawindow,1);
numcells = size(fields,1);
windowcount = 0;

for k = 1:numwindows
    spikes = cell(numcells,2);
    spkcount = 0;
    for i = 1:numcells
        spkDS = [];
        spktimes = fields {i,5}; %get all spiketimes for cell i
        spkDS = spktimes(spktimes>=thetawindow(k,1) & spktimes<=thetawindow(k,2));%get the ones that fall in between the window
        for j = 1:size(fields{i,2},2)
           spkDS(:,j+1) = fields{i,2}(1,j); 
        end
        if size(spkDS)>0
        spikes{i,1} = spkDS(:,1);
        spikes{i,2} = spkDS(:,2:end);
        else
        spikes{i,1} = [];
        spikes{i,2} = [];
        end
        if size(spkDS,1) > 0
        spkcount = spkcount + 1; %for limiting number of active cells
        %spkcount = spkcount + size(spkDS(:,1),1); %for limiting total spikes
        end
    end
    if spkcount >= minspikes
        windowcount = windowcount +1;
        result{windowcount,1} = spikes;
        result{windowcount,2} = thetawindow(k,:);
    end
end

% for i = 1:numcells
%    if ~isempty(sp{i,1}) 
%    times(:,1) = [times(:,1) ; sp{i,1}];
%    timesL(:,1) = [timesL(:,1) ; sp{i,2}(:,1)];
%    end
% end