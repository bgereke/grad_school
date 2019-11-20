function [map, xAxis] = probX(posx, binWidth)

%this version includes the smoothing from raymonds episode_1_8

% Number of bins of the map
tLength = max(posx) - min(posx);
numBins = ceil(tLength / binWidth);

% Allocate memory for the maps
timeMap = zeros((numBins),1);
xAxis = zeros((numBins),1);

startPos =  min(posx);
stopPos = startPos + binWidth;
 
for ii = 1:1:numBins
    timeMap(ii) = length(find(posx >= startPos & posx < stopPos));
    % Set the axis
    xAxis(ii) = startPos+ binWidth/2;
        
    % Increment the position bin
    startPos = startPos + binWidth;
    stopPos = stopPos + binWidth;
end
   map = timeMap./sum(timeMap);


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
