% Load video data and interpolate missing samples
% Example of how to load data
%[x,y,t] = readVideoData(videoFile, positionScaleFactor);

function [x,y,t,a] = readVideoData(file,scale)


% Set the field selection for reading the video files. 1 = Add parameter, 
% 0 = skip parameter
fieldSelection(1) = 1; % Timestamps
fieldSelection(2) = 1; % Extracted X
fieldSelection(3) = 1; % Extracted Y
fieldSelection(4) = 1; % Extracted Angel
fieldSelection(5) = 0; % Targets
fieldSelection(6) = 0; % Points
% Do we return header 1 = Yes, 0 = No.
extractHeader = 0;
% 5 different extraction modes, see help file for Nlx2MatVT
extractMode = 1; % Extract all data

% Use the NeuraLynx mex-compiled file to read the video data into the
% Matlab work memory.
[t, x, y, a] = Nlx2MatVT(file,fieldSelection,extractHeader,extractMode);

% Convert timestamps to seconds
t = t/1000000;

% Interpolate missing segments of the path, but max gap of 4 seconds
[x,y] = interporPos(x,y,4,25);

% Set the (0,0)-coordinates to NaN because these are missing position
% samples
ind = find(x==0 & y==0);
x(ind) = NaN;
y(ind) = NaN;
x = inpaint_nans(x);
y = inpaint_nans(y);

% Transform the coordinates from pixels to centimeters
x = x * scale;
y = y * scale;

% Smooth the path with a moving mean window filter (removing tracker
% jitter)
for cc=8:length(x)-7
    x(cc) = nanmean(x(cc-7:cc+7));
    y(cc) = nanmean(y(cc-7:cc+7));
end

