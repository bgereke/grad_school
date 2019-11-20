function [times] = timeQuadrants3(startMinute,endMinute)

%this version doesn't use a "corners" file, instead taking the max and min
%of the positions to establish the environment, so make sure the coverage
%was good before using.

%need to figure out scale still
scaler = 1;

%Positions for current folder which should be beginx
VTfile = ('./VT1.nvt');
[timepoints,xpoints,ypoints] = getPositions(VTfile,scaler);

%truncate data for desired time
ind = find (  (timepoints >= ( min(min(timepoints)) + startMinute*60)) ...
            & (timepoints < ( min(min(timepoints)) + endMinute*60  )) );
    
posT = timepoints(ind);
posX = xpoints (ind);
posY = ypoints (ind);

%get corners of box for the day using the whole sessions position file
xmax = max(max(xpoints));
xmin = min(min(xpoints));
ymax = max(max(ypoints));
ymin = min(min(ypoints));
xspan = xmax - xmin;
yspan = ymax - ymin;
vlineX = [(xmin + xspan/2) (xmin + xspan/2)];
vlineY = [(ymax) (ymin)];
hlineX = [(xmin) (xmax)];
hlineY = [(ymin + yspan/2) (ymin + yspan/2)];
posY = -posY + yspan + 2*ymin; % reflects posisitons so they are consistent with orientation in recording room
ypoints = -ypoints + yspan + 2*ymin;

%plot 
figure;
hold on;
plot (xpoints,ypoints); 
plot (posX,posY,'r','linewidth',2);
plot (xmin, ymin,'o', xmin, ymax,'o', xmax, ymin,'o', xmax, ymax,'o');
plot (vlineX,vlineY, 'linewidth',2);
plot (hlineX,hlineY, 'linewidth',2);

%tally up time points from original file that land within the four main
%bins

totDur = posT(end) - posT(1);
sampDur = totDur/length(posX);

times = zeros(4);

ind = find (posX < (xmin+xspan/2) ...
          & posY > (ymin+yspan/2));
times (1,1) = length(ind)*(sampDur);

ind = find (posX > (xmin+xspan/2) ...
          & posY > (ymin+yspan/2));
times (1,2) = length(ind)*(sampDur);

ind = find (posX < (xmin+xspan/2) ...
          & posY < (ymin+yspan/2));
times (2,1) = length(ind)*(sampDur);

ind = find (posX > (xmin+xspan/2) ...
          & posY < (ymin+yspan/2));
times (2,2) = length(ind)*(sampDur);





function [t,x,y] = getPositions(file,scale)
% Set the field selection for reading the video files. 1 = Add parameter, 
% 0 = skip parameter
fieldSelection(1) = 1; % Timestamps
fieldSelection(2) = 1; % Extracted X
fieldSelection(3) = 1; % Extracted Y
fieldSelection(4) = 0; % Extracted Angel
fieldSelection(5) = 0; % Targets
fieldSelection(6) = 0; % Points
% Do we return header 1 = Yes, 0 = No.
extractHeader = 0;
% 5 different extraction modes, see help file for Nlx2MatVT
extractMode = 1; % Extract all data

% Use the NeuraLynx mex-compiled file to read the video data into the
% Matlab work memory.
[t, x, y] = Nlx2MatVT(file,fieldSelection,extractHeader,extractMode);

% Convert timestamps to seconds
t = t/1000000;

% Interpolate missing segments of the path, but max gap of 4 seconds
[x,y] = interporPos(x,y,4,25);

% Set the (0,0)-coordinates to NaN because these are missing position
% samples
ind = find(x==0 & y==0);
x(ind) = NaN;
y(ind) = NaN;

% Transform the coordinates from pixels to centimeters
x = x * scale;
y = y * scale;

%centers positions to 0,0
%[x,y] = centreBox(x,y);

% Smooth the path with a moving mean window filter (removing tracker
% jitter)
for cc=8:length(x)-7
    x(cc) = nanmean(x(cc-7:cc+7));
    y(cc) = nanmean(y(cc-7:cc+7));
end


