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