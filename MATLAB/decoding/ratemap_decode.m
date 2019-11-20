function [map] = ratemap_decode(posx, spkx, vel,spkvel, wsgzmu, wfgzmu, spksfgz, win, xbins)

% Allocate memory for the maps
numBins = length(xbins);
spikeMap = zeros((numBins),2);
timeMap = zeros((numBins),2);
 
for ii = 1:1:numBins
    spikeMap(ii,1) = length(spkx(spkvel>=0 & spkx == xbins(ii)));    
    spikeMap(ii,2) = length(spkx(spkvel<0 & spkx == xbins(ii)));    
    spikeMap(ii,3) = length(spkx(spkvel>=0 & spkx == xbins(ii) & spksfgz(:,1)>0.25 & spksfgz(:,2)<0.5));
    spikeMap(ii,4) = length(spkx(spkvel<0 & spkx == xbins(ii) & spksfgz(:,1)>0.25 & spksfgz(:,2)<0.5));
    spikeMap(ii,5) = length(spkx(spkvel>=0 & spkx == xbins(ii) & spksfgz(:,1)<0.5 & spksfgz(:,2)>0.25));
    spikeMap(ii,6) = length(spkx(spkvel<0 & spkx == xbins(ii) & spksfgz(:,1)<0.5 & spksfgz(:,2)>0.25)); 
    
    timeMap(ii,1) = sum(win(posx==xbins(ii) & vel>=0,2)-win(posx==xbins(ii) & vel>=0,1));
    timeMap(ii,2) = sum(win(posx==xbins(ii) & vel<0,2)-win(posx==xbins(ii) & vel<0,1));
    timeMap(ii,3) = sum(win(posx==xbins(ii) & vel>=0 & wsgzmu>0.25 & wfgzmu<0.5,2)...
        -win(posx==xbins(ii) & vel>=0 & wsgzmu>0.25 & wfgzmu<0.5,1));
    timeMap(ii,4) = sum(win(posx==xbins(ii) & vel<0 & wsgzmu>0.25 & wfgzmu<0.5,2)...
        -win(posx==xbins(ii) & vel<0 & wsgzmu>0.25 & wfgzmu<0.5,1));
    timeMap(ii,5) = sum(win(posx==xbins(ii) & vel>=0 & wsgzmu<0.5 & wfgzmu>0.25,2)...
        -win(posx==xbins(ii) & vel>=0 & wsgzmu<0.5 & wfgzmu>0.25,1));
    timeMap(ii,6) = sum(win(posx==xbins(ii) & vel<0 & wsgzmu<0.5 & wfgzmu>0.25,2)...
        -win(posx==xbins(ii) & vel<0 & wsgzmu<0.5 & wfgzmu>0.25,1));
end

    % Smooth the spike and time maps
    spikeMap(:,1) = mapSmoothing(spikeMap(:,1));
    timeMap(:,1) = mapSmoothing(timeMap(:,1));
    spikeMap(:,2) = mapSmoothing(spikeMap(:,2));    
    timeMap(:,2) = mapSmoothing(timeMap(:,2));
    spikeMap(:,3) = mapSmoothing(spikeMap(:,3));    
    timeMap(:,3) = mapSmoothing(timeMap(:,3));
    spikeMap(:,4) = mapSmoothing(spikeMap(:,4));    
    timeMap(:,4) = mapSmoothing(timeMap(:,4));
    spikeMap(:,5) = mapSmoothing(spikeMap(:,5));    
    timeMap(:,5) = mapSmoothing(timeMap(:,5));
    spikeMap(:,6) = mapSmoothing(spikeMap(:,6));    
    timeMap(:,6) = mapSmoothing(timeMap(:,6));

    timeMap(timeMap<0.5) = 0;
    %timeMap(timeMap < 1) = 0;
    map = spikeMap ./ timeMap;


% Smooths the map with guassian smoothing
function sMap = mapSmoothing(map)

box = boxcarTemplate1D();

numBins = length(map);
sMap = zeros(numBins,1);

for ii = 1:numBins
    for k = 1:5
        % Bin shift
        sii = k-3;
        % Bin index
        binInd = ii + sii;
        % Boundry check
        if binInd<1
            binInd = 1;
        end
        if binInd > numBins
            binInd = numBins;
        end
        
        sMap(ii) = sMap(ii) + map(binInd) * box(k);
    end
end

% 1-D Gaussian boxcar template
function box = boxcarTemplate1D()

% Gaussian boxcar template
box = [0.05 0.25 0.40 0.25 0.05];
