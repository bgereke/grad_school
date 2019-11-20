function [Pvx,Pxv,Pv,Px] = probXV(Px,posx,v,windowV,binWidth)

%calculate probability of position x, given velocity v (Pxv)
%need Px (provided), Pv, and Pvx
%get window velocities
vmin = windowV - 1000;%set really high to simply control for direction
vmax = windowV + 1000;
if windowV < 0 && vmax > 0 %comment out if do NOT want Pvx to reveal direction, when the previous is set really high
    vmax = 0;
elseif windowV > 0 && vmin < 0
    vmin = 0;
end

%probabilty of v within +/- 3 cm/sec
Pv = length(find(v>vmin & v<vmax)) / length(v);

%find Pvx at each x bin
% Number of bins of the map
tLength = max(posx) - min(posx);
numBins = ceil(tLength / binWidth);

% Allocate memory for the maps
startPos =  min(posx);
stopPos = startPos + binWidth;
 
for ii = 1:1:numBins
    indx = find(posx >= startPos & posx < stopPos);
    tempv = v(indx);
    Pvx(ii,1) = length(find(tempv>vmin & tempv<vmax)) / length(tempv);
        
    % Increment the position bin
    startPos = startPos + binWidth;
    stopPos = stopPos + binWidth;
end

Pxv = Pvx.*Px./Pv;




