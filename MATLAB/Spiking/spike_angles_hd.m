% equalPlot_max('inFile.txt',scale_max)
%
% equalPlot generates ratemaps, plots them and store them to images.
% Multiple session will be read and the maximum rate for a cell will be
% used for plotting in all sessions.
%
% The input file must be on the following format.
%
% C:\Data\TTList.txt
% C:\Data\Begin 1
% C:\Data\Begin 2
% C:\Data\Begin 3
% C:\Data\Begin 4
% and so on ...
%
% The first line specifie the directory to the first session and also the
% name on the t-file list that will be used for all the sessions listed.
% The t-file list must contain the Mclust t-files, one name for each file.
% If the list contain cells that only occure in some of the sessions, these
% cells will be plotted as having zero firing rate when missing. Placefield
% images will be stored to both bmp and eps imagefiles to a subdirectory in
% the data folder called placeFieldImages.

function [peakRate, pathCoord] = spike_angles_hd(inFile,scale_max)

D = 100; % Circle diameter in cm
%bins = 60; % Number of bins (roughly 3.5 cm)
%scale = 0.26;% Conversion from pixels to cm, frame is 640 by 480 (?)
%scale = 0.5;% Conversion from pixels to cm, frame is 640 by 480 (?)
img_text = 'on';

fid = fopen(inFile,'r');
if fid == -1
    msgbox('Could not open the input file! Make sure the filname and path are correct.','ERROR');
end

% Get sessions and csc-file list from input file
fid = fopen(inFile,'r');
ii = -1;
numsessions = 0;
while ~feof(fid)
    str = fgetl(fid);
    if ii == -1
        ttList = str;
    elseif ii == 0
        cscList = str;
    elseif ii == 1
        refList = str;
    elseif ii == 2
        ch4avg = str;
    elseif ii > 2
        numsessions  = numsessions+1;
        if ~strcmp(str(end),'\')
            str = strcat(str,'\');
        end
        sessions(numsessions) = {str};
    end
    ii = ii+1;
end
counter = numsessions;

% read the file names from the t-file list
F = ReadFileList(ttList);
% Number of cells/files in the list
numCells = length(F);

% Set the field selection for reading the video files. 1 = Add parameter, 0 = skip
% parameter
fieldSelection(1) = 1; % Timestamps
fieldSelection(2) = 0; % Extracted X
fieldSelection(3) = 0; % Extracted Y
fieldSelection(4) = 0; % Extracted Angel
fieldSelection(5) = 1; % Targets
fieldSelection(6) = 0; % Points
% Do we return header 1 = Yes, 0 = No.
extractHeader = 0;
% 5 different extraction modes, see help file for Nlx2MatVt
extractMode = 1; % Extract all data

% Read the input data
%binWidth = pi*D/bins;
pathCoord = cell(counter,4);
spk_phase = cell(counter,numCells);
ts = cell(counter,numCells);
map = cell(counter,numCells);
phase = cell(counter,1);
t = cell(counter,1);
hd = cell(counter,1);
vel = cell(counter,1);

for ii = 1:counter
    disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));
    % Check if subdir for storing images are present. If not, it is
    % created.
    temp = sessions{ii};
    dirInfo = dir(temp);
    found1 = 0;
    found2 = 0;
    for kk=1:size(dirInfo,1)
        if dirInfo(kk).isdir
            if strcmp(dirInfo(kk).name,strcat('SpikeAngles_hd','\'))
                found1 = 1;
            end
            if strcmp(dirInfo(kk).name,strcat('PlaceFieldAngles_hd','\'))
                found2 = 1;
            end
        end
    end
    if found1==0
        mkdir(temp,strcat('SpikeAngles_hd','\'));
    end
    if found2==0
        mkdir(temp,strcat('PlaceFieldAngles_hd','\'));
    end
    
    % Get position and head direction data
    file = strcat(sessions{ii},'vt1.nvt');
    [t{ii}, targets] = Nlx2MatVT(file,fieldSelection,extractHeader,extractMode);
    % Decode the target data
    [dTargets,tracking] = decodeTargets(targets);
    % Exctract position data from the target data
    [frontX,frontY,backX,backY] = extractPosition(dTargets,tracking);
    % Interpolate missing position samples
    [frontX] = interporPos(frontX,1.5,29.97);[frontY] = interporPos(frontY,1.5,29.97);
    [backX] = interporPos(backX,1.5,29.97);[backY] = interporPos(backY,1.5,29.97);
    % Smooth the exctracted position samples using a moving mean filter
    % Set missing samples to NaN. Necessary when using the mean filter.
    indf = find(frontX==0);
    frontX(indf) = NaN;
    frontY(indf) = NaN;
    indb = find(backX==0);
    backX(indb) = NaN;
    backY(indb) = NaN;
    [frontX,frontY] = posMeanFilter(frontX,frontY);
    [backX,backY] = posMeanFilter(backX,backY);
    x = nanmean([frontX;backX]);
    y = nanmean([frontY;backY]);
    % Set the missing samples back to zero
    frontX(indf) = 0;x(indf) = 0;
    frontY(indf) = 0;y(indf) = 0;
    backX(indb) = 0;x(indb) = 0;
    backY(indb) = 0;y(indb) = 0;
    
    ind = find(x == 0);
    t{ii}(ind) = [];
    x(ind) = [];frontX(ind) = [];backX(ind) = [];
    y(ind) = [];frontY(ind) = [];backY(ind) = [];
    
    %find track edges
    [hx,xc] = hist(x,50); [hy,yc] = hist(y,50);
    th = 50;
    xmin = xc(find(hx>th,1,'first'));
    xmax = xc(find(hx>th,1,'last'));
    ymin = yc(find(hy>th,1,'first'));
    ymax = yc(find(hy>th,1,'last'));
    
    % Adjust input data
    scalex = D/(xmax-xmin);
    scaley = D/(ymax-ymin);
    x = x * scalex;frontX = frontX*scalex;backX = backX*scalex;
    y = y * scaley;frontY = frontY*scaley;backY = backY*scaley;
    xmin = xmin * scalex;
    ymin = ymin * scaley;
    xmax = xmax * scalex;
    ymax = ymax * scaley;
    xcen = xmin+(xmax-xmin)/2;
    ycen = ymin+(ymax-ymin)/2;
    x = x - xcen;frontX = frontX-xcen;backX = backX-xcen;
    y = y - ycen;frontY = frontY-ycen;backY = backY-ycen;
    phase{ii} = atan2(y,x);
    % Calculate the head stage direction (-pi to pi)
    hd{ii} = headDirection(frontX,frontY,backX,backY);
    hd{ii} = hd{ii}/360*2*pi - pi;
    % determine direction relitive to track
    rot = [cos(-pi/2),-sin(-pi/2);sin(-pi/2),cos(-pi/2)];
    track = (rot*[cos(phase{ii})',sin(phase{ii})']')';
    track = atan2(track(:,2),track(:,1));
    hd{ii} = angle(exp(1j*track).*conj(exp(1j*hd{ii})))';
    % rotate if red light wasn't in front
    if abs(circ_mean(hd{ii}')) > pi/2
        hd{ii} = (rot*[cos(hd{ii})',sin(hd{ii})']')';
        hd{ii} = atan2(hd{ii}(:,2),hd{ii}(:,1));
    end
    % Convert timestamps to seconds
    t{ii} = t{ii}/1000000;
    
    %compute velocity and velocity map
    vel{ii} = zeros(1,length(x));
    xs = smooth(x,15);ys = smooth(y,15);
    for v = 2:1:(length(x)-1)
        vel{ii}(v) = sqrt((xs(v+1)-xs(v-1))^2 + (ys(v+1)-ys(v-1))^2)/(t{ii}(v+1)-t{ii}(v-1));
    end
    vel{ii}(end) = sqrt((xs(end)-xs(end-1))^2 + (ys(end)-ys(end-1))^2)/(t{ii}(length(x))-t{ii}(length(x)-1));
    vel{ii}(vel{ii}>=40) = 0.5*(vel{ii}(circshift((vel{ii}>=40),-3))+vel{ii}(circshift((vel{ii}>=40),3)));
    vel{ii}(:) = smooth(vel{ii},15);
    
    % Store path
    pathCoord{ii,1} = x;
    pathCoord{ii,2} = y;
    pathCoord{ii,3} = t{ii};
    pathCoord{ii,4} = phase{ii};
    
    % END COMPUTING POSITION DATA
    
    % Open TTL Event Data
    efile = strcat(sessions{ii},'Events.nev');
    FS = [1 0 1 0 1]; EH = 0; EM = 1;
    [TimeStamps{ii}, TTLs{ii}, EventStrings] = Nlx2MatEV(efile,FS,EH,EM);
    TTLs{ii}(1) = []; TimeStamps{ii}(1) = []; TimeStamps{ii} = TimeStamps{ii}/1000000;
    
    % Read the spike data, but must look out for files in the list that
    % don't exist for this session
    dirFiles = dir(sessions{ii});   % Files in directory
    NaFile = cell(1,length(F)); % Store non existing file names here
    NaFileCounter = 0;          % Counts number of non existing files
    for cc=1:length(F)
        I = strmatch(char(F(cc)),char(dirFiles.name),'exact'); % Try to find file in dir
        if length(I)==0 % File is non existing
            NaFileCounter = NaFileCounter + 1;
            NaFile(1,NaFileCounter) = F(cc);
        end
    end
    NaFile = NaFile(1,1:NaFileCounter);
    
    binsize = 5; %cm
    fieldAxis = linspace(-pi,pi,pi*D/binsize);
    h = 100; % Smoothing factor when calculating the ratemap
    minr = 1; minp = 0.2; %min rate and proportion of peak rate
    
    % Load data from the cut files generated by Mclust
    S = LoadSpikes(F,sessions{ii},NaFile);
    disp('Start calculating the ratemaps for the cells');
    for jj=1:numCells
        disp(sprintf('%s%i',' Cell ',jj));
        if ~isa(S{jj},'ts') % Empty cell in this session
            map{ii,jj} = zeros(floor(pi*D/binsize));
            ts{ii,jj} = 1e64; % use a single ridicilous time stamp if the cell is silent
        else
            % Convert t-file data to timestamps in second
            ts{ii,jj} = Data(S{jj}) / 10000;
            % Get position to spikes
            [spk_phase{ii,jj}] = spikePos(ts{ii,jj},phase{ii},t{ii});
            [spk_vel{ii,jj}] = spikePos(ts{ii,jj},vel{ii},t{ii});
        end
    end
end
disp('Plot maps and store them to files');

for ii = 1:counter
    TimeStamps{ii} = TimeStamps{ii} - min(t{ii});
    t{ii} = t{ii}-min(t{ii});
    for jj=1:numCells
        close all
        figure(1); hold on
        if length(TTLs{ii})>=2
            if TTLs{ii}(1) == 1
                rectangle('Position',[min(t{ii}) min(phase{ii}),TimeStamps{ii}(1),max(phase{ii})-min(phase{ii})],'FaceColor','g','EdgeColor','g')
            end
            for i = 2:length(TTLs{ii})
                if TTLs{ii}(i) == 1
                    rectangle('Position',[TimeStamps{ii}(i-1),min(phase{ii}),TimeStamps{ii}(i)-TimeStamps{ii}(i-1),max(phase{ii})-min(phase{ii})],'FaceColor','g','EdgeColor','g')
                end
                tidx = length(t{ii}(t{ii}<=TimeStamps{ii}(i)))+1;
            end
            if TTLs{ii}(end) == 0
                rectangle('Position',[TimeStamps{ii}(end),min(phase{ii}),max(t{ii})-TimeStamps{ii}(end),max(phase{ii})-min(phase{ii})],'FaceColor','g','EdgeColor','g')
            end
        end
        z = zeros(size(t{ii}));
        surface([t{ii};t{ii}],[phase{ii};phase{ii}],[z;z],[hd{ii};hd{ii}],'facecol','no','edgecol','interp','linew',6); %plot trajectory w/hd specifying color
        idx = zeros(length(spk_phase{ii,jj}),1);
        for s = 1:length(spk_phase{ii,jj})
            [~,idx(s)] = min((phase{ii}-spk_phase{ii,jj}(s)).^2);
        end
        plot(t{ii}(idx),phase{ii}(idx),'.k');title(strcat('Cell ',num2str(jj),'-',F{jj}))
        xlabel('time (sec)');ylabel('track angle');colorbar;colormap jet
        caxis([-pi/2 pi/2])
        
        bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('SpikeAngles_hd','\'),F{jj}(1:end-2),'.bmp');
        epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('SpikeAngles_hd','\'),F{jj}(1:end-2),'.eps');
        f = getframe(gcf);
        [pic, cmap] = frame2im(f);
        imwrite(pic,bmpImage,'bmp');
        saveas(1,epsImage,'epsc');
        pause(0.1)
        
    end
end



%__________________________________________________________________________
%
% Field functions
%__________________________________________________________________________

% timeMap will find the amount of time the rat spends in each bin of the
% box. The box is divided into the same number of bins as in the rate map.
% posx, posy, post: Path coordinates and time stamps.
% bins: Number of bins the box is diveded int (bins x bins)
% extremal: minimum and maximum positions.
function timeMap = findTimeMap(phase,post,bins,binWidth)

% Duration of trial
duration = post(end)-post(1);
% Average duration of each position sample
sampDur = duration/length(phase);

% phase for current bin position
pcx = -pi-binWidth;
timeMap = zeros(bins);
% Find number of position samples in each bin
for i = 1:bins
    % Increment the pahse
    pcx = pcx + binWidth;
    I = find(phase >= pcx & phase < pcx+binWidth);
    % Number of position samples in the current bin
    timeMap(i) = length(I);
end
% Convert to time spent in bin
timeMap = timeMap * sampDur;


% Calculates the rate map.
function map = ratemap(ts,spk_phase,phase,post,h,mapAxis)
invh = 1/h;
map = zeros(length(mapAxis));
idx = 0;
for i = mapAxis
    idx = idx + 1;
    map(idx) = rate_estimator(ts,spk_phase,i,invh,phase,post);
end

% Calculate the rate for one position value
function r = rate_estimator(ts,spk_phase,i,invh,phase,post)
% edge-corrected kernel density estimator
% make (spk_phase-i)->-pi:pi
delta = min([abs(spk_phase - i),2*pi-abs(spk_phase-i)],[],2);
conv_sum = sum(gaussian_kernel(delta*invh));
delta = min([abs(phase - i);2*pi-abs(phase-i)],[],1);
edge_corrector =  trapz(post,gaussian_kernel((delta*invh)));
%edge_corrector(edge_corrector<0.15) = NaN;
r = (conv_sum / (edge_corrector + 0.1)) + 0.1; % regularised firing rate for "wellbehavedness"
% i.e. no division by zero or log of zero
% Gaussian kernel for the rate calculation
function r = gaussian_kernel(x)
% k(u) = ((2*pi)^(-length(u)/2)) * exp(u'*u)
r = 0.15915494309190 * exp(-0.5*(x.*x));


% Finds the position to the spikes
function [spk_phase] = spikePos(ts,phase,post)
N = length(ts);
spk_phase = zeros(N,1);
for ii = 1:N
    tdiff = (post-ts(ii)).^2;
    [~,ind] = min(tdiff);
    spk_phase(ii) = phase(ind(1));
end

%__________________________________________________________________________
%
% Function for modifing the path
%__________________________________________________________________________

% Median filter for positions
function [posx,posy,post] = medianFilter(posx,posy,post)
N = length(posx);
x = [NaN*ones(15,1); posx'; NaN*ones(15,1)];
y = [NaN*ones(15,1); posy'; NaN*ones(15,1)];
X = ones(1,N);
Y = ones(1,N);
for cc=16:N+15
    lower = cc-15;
    upper = cc+15;
    X(cc-15) = nanmedian(x(lower:upper));
    Y(cc-15) = nanmedian(y(lower:upper));
end
index = find(isfinite(X));
posx=X(index);
posy=Y(index);
post=post(index);

% Removes position values that are a result of bad tracking
function [posx,posy,post] = offValueFilter(posx,posy,post)

x1 = posx(1:end-1);
x2 = posx(2:end);
y = posy(2:end);
t = post(2:end);

N = 1;
while N > 0
    len = length(x1);
    dist = abs(x2 - x1);
    % Finds were the distance between two neighbour samples are more than 2
    % cm.
    ind = find(dist > 2);
    N = length(ind);
    if N > 0
        if ind(1) ~= len
            x2(:,ind(1)) = [];
            x1(:,ind(1)+1) = [];
            y(:,ind(1)) = [];
            t(:,ind(1)) = [];
        else
            x2(:,ind(1)) = [];
            x1(:,ind(1)) = [];
            y(:,ind(1)) = [];
            t(:,ind(1)) = [];
        end
    end
end
x = x2(2:end);
t = t(2:end);
y1 = y(1:end-1);
y2 = y(2:end);

N = 1;
while N > 0
    len2 = length(y1);
    dist = abs(y2 - y1);
    ind = find(dist > 3);
    N = length(ind);
    if N > 0
        if ind(1) ~= len2
            y2(:,ind(1)) = [];
            y1(:,ind(1)+1) = [];
            x(:,ind(1)) = [];
            t(:,ind(1)) = [];
        else
            y2(:,ind(1)) = [];
            y1(:,ind(1)) = [];
            x(:,ind(1)) = [];
            t(:,ind(1)) = [];
        end
    end
end
posx = x;
posy = y2;
post = t;

% Centre the path/box for the two sessions. Both boxes are set according to
% the path/coordinates to the box for the first session.
function [posx,posy] = centreBox(posx,posy)

% Find border values for box for session 1
maxX = max(posx);
minX = min(posx);
maxY = max(posy);
minY = min(posy);

% Set the corners of the reference box
NE = [maxX, maxY];
NW = [minX, maxY];
SW = [minX, minY];
SE = [maxX, minY];

% Get the centre coordinates of the box
centre = findCentre(NE,NW,SW,SE);
% Centre both boxes according to the coordinates to the first box
posx = posx - centre(1);
posy = posy - centre(2);


% Calculates the centre of the box from the corner coordinates
function centre = findCentre(NE,NW,SW,SE);

% The centre will be at the point of interception by the corner diagonals
a = (NE(2)-SW(2))/(NE(1)-SW(1)); % Slope for the NE-SW diagonal
b = (SE(2)-NW(2))/(SE(1)-NW(1)); % Slope for the SE-NW diagonal
c = SW(2);
d = NW(2);
x = (d-c+a*SW(1)-b*NW(1))/(a-b); % X-coord of centre
y = a*(x-SW(1))+c; % Y-coord of centre
centre = [x,y];


%__________________________________________________________________________
%
% Additional graphics function
%__________________________________________________________________________

function drawfield(map,faxis,cmap,maxrate,cellid,cell_file,text)

% This function will calculate an RGB image from the rate
% map. We do not just call image(map) and caxis([0 maxrate]),
% as it would plot unvisted parts with the same colour code
% as 0 Hz firing rate. Instead we give unvisited bins
% their own colour (e.g. gray or white).

maxrate = ceil(maxrate);
if maxrate < 1
    maxrate = 1;
end
n = size(map,1);
plotmap = ones(n,n,3);
for jj = 1:n
    for ii = 1:n
        if isnan(map(jj,ii))
            plotmap(jj,ii,1) = 1; % give the unvisited bins a gray colour
            plotmap(jj,ii,2) = 1; %KJ changed from 1 1 1
            plotmap(jj,ii,3) = 1;
        else
            if (map(jj,ii) > maxrate)
                plotmap(jj,ii,1) = 1; % give the unvisited bins a gray colour
                plotmap(jj,ii,2) = 0;
                plotmap(jj,ii,3) = 0;
            else
                rgb = pixelcolour(map(jj,ii),maxrate,cmap);
                plotmap(jj,ii,1) = rgb(1);
                plotmap(jj,ii,2) = rgb(2);
                plotmap(jj,ii,3) = rgb(3);
            end
        end
    end
end
image(faxis,faxis,plotmap);
set(gca,'YDir','Normal');
axis image;
axis off
if strcmp(text,'on')
    title(strcat(cellid,'¤(0 - ',num2str(maxrate), ' Hz)','¤',cell_file),'FontSize',20);
end


function rgb = pixelcolour(map,maxrate,cmap)

% This function calculates a colour for each bin
% in the rate map.

cmap1 = ...
    [    0         0    0.5625; ...
    0         0    0.6875; ...
    0         0    0.8125; ...
    0         0    0.9375; ...
    0    0.0625    1.0000; ...
    0    0.1875    1.0000; ...
    0    0.3125    1.0000; ...
    0    0.4375    1.0000; ...
    0    0.5625    1.0000; ...
    0    0.6875    1.0000; ...
    0    0.8125    1.0000; ...
    0    0.9375    1.0000; ...
    0.0625    1.0000    1.0000; ...
    0.1875    1.0000    0.8750; ...
    0.3125    1.0000    0.7500; ...
    0.4375    1.0000    0.6250; ...
    0.5625    1.0000    0.5000; ...
    0.6875    1.0000    0.3750; ...
    0.8125    1.0000    0.2500; ...
    0.9375    1.0000    0.1250; ...
    1.0000    1.0000         0; ...
    1.0000    0.8750         0; ...
    1.0000    0.7500         0; ...
    1.0000    0.6250         0; ...
    1.0000    0.5000         0; ...
    1.0000    0.3750         0; ...
    1.0000    0.2500         0; ...
    1.0000    0.1250         0; ...
    1.0000         0         0; ...
    0.8750         0         0; ...
    0.7500         0         0; ...
    0.6250         0         0 ];

cmap2 = ...
    [0.0417         0         0; ...
    0.1250         0         0; ...
    0.2083         0         0; ...
    0.2917         0         0; ...
    0.3750         0         0; ...
    0.4583         0         0; ...
    0.5417         0         0; ...
    0.6250         0         0; ...
    0.7083         0         0; ...
    0.7917         0         0; ...
    0.8750         0         0; ...
    0.9583         0         0; ...
    1.0000    0.0417         0; ...
    1.0000    0.1250         0; ...
    1.0000    0.2083         0; ...
    1.0000    0.2917         0; ...
    1.0000    0.3750         0; ...
    1.0000    0.4583         0; ...
    1.0000    0.5417         0; ...
    1.0000    0.6250         0; ...
    1.0000    0.7083         0; ...
    1.0000    0.7917         0; ...
    1.0000    0.8750         0; ...
    1.0000    0.9583         0; ...
    1.0000    1.0000    0.0625; ...
    1.0000    1.0000    0.1875; ...
    1.0000    1.0000    0.3125; ...
    1.0000    1.0000    0.4375; ...
    1.0000    1.0000    0.5625; ...
    1.0000    1.0000    0.6875; ...
    1.0000    1.0000    0.8125; ...
    1.0000    1.0000    0.9375];

if strcmp(cmap,'jet')
    steps = (31*(map/maxrate))+1;
    steps = round(steps);
    rgb = cmap1(steps,:);
else
    steps = (31*(map/maxrate))+1;
    steps = round(steps);
    rgb = cmap2(steps,:);
end

%__________________________________________________________________________
%
%      Functions for reading Mclust data
%__________________________________________________________________________

function S = LoadSpikes(tfilelist, path, NaFile)
% tfilelist:    List of t-files. Each file contains a cluster of spikes
%               from a cell.
% path:         Path to the directory were the t-files are stored
% NaFile:       List of file names in tfilelist that don't exist
%               in the current directory
%
% inp: tfilelist is a cellarray of strings, each of which is a
% tfile to open.  Note: this is incompatible with version unix3.1.
% out: Returns a cell array such that each cell contains a ts
% object (timestamps which correspond to times at which the cell fired)
%
% Edited by: Raymond Skjerpeng


%-------------------
% Check input type
%-------------------

if ~isa(tfilelist, 'cell')
    error('LoadSpikes: tfilelist should be a cell array.');
end


% Number of file names in tfilelist
nFiles = length(tfilelist);
% Actual number of files to be loaded
anFiles = nFiles - length(NaFile);


%--------------------
% Read files
%--------------------
fprintf(2, 'Reading %d files.', anFiles);

% for each tfile
% first read the header, then read a tfile
% note: uses the bigendian modifier to ensure correct read format.


S = cell(nFiles, 1);
for iF = 1:nFiles
    %DisplayProgress(iF, nFiles, 'Title', 'LoadSpikes');
    tfn = tfilelist{iF};
    % Check if file exist
    if length(strmatch(char(tfn),NaFile,'exact'))>0
        S{iF} = -1; % Set this as default for each file that doesn't exist
    else
        tfn = strcat(strcat(path,'\'),tfn); % Path to file + file name
        if ~isempty(tfn)
            tfp = fopen(tfn, 'rb','b');
            if (tfp == -1)
                warning([ 'Could not open tfile ' tfn]);
            end
            
            ReadHeader(tfp);
            S{iF} = fread(tfp,inf,'uint32');	%read as 32 bit ints
            S{iF} = ts(S{iF});
            fclose(tfp);
        end 	% if tfn valid
    end
end		% for all files
fprintf(2,'\n');


function F = ReadFileList(fn)

% F = ReadFileList(fn)
%
% INPUTS:
%   fn -- an ascii file of filenames, 1 filename per line
%
% OUTPUTS:
%   F -- a cell array of filenames suitable for use in programs
%        such as LoadSpikes
%
% Now can handle files with headers

% ADR 1998
% version L4.1
% status: PROMOTED

% v4.1 added ReadHeader

[fp,errmsg] = fopen(fn, 'rt');
if (fp == -1)
    error(['Could not open "', fn, '". ', errmsg]);
end

ReadHeader(fp);
ifp = 1;
while (~feof(fp))
    F{ifp} = fgetl(fp);
    ifp = ifp+1;
end
fclose(fp);

F = F';


function H = ReadHeader(fp)
% H = ReadHeader(fp)
%  Reads NSMA header, leaves file-read-location at end of header
%  INPUT:

%      fid -- file-pointer (i.e. not filename)
%  OUTPUT:
%      H -- cell array.  Each entry is one line from the NSMA header
% Now works for files with no header.
% ADR 1997
% version L4.1
% status: PROMOTED
% v4.1 17 nov 98 now works for files sans header
%---------------

% Get keys
beginheader = '%%BEGINHEADER';
endheader = '%%ENDHEADER';

iH = 1; H = {};
curfpos = ftell(fp);

% look for beginheader
headerLine = fgetl(fp);
if strcmp(headerLine, beginheader)
    H{1} = headerLine;
    while ~feof(fp) & ~strcmp(headerLine, endheader)
        headerLine = fgetl(fp);
        iH = iH+1;
        H{iH} = headerLine;
    end
else % no header
    fseek(fp, curfpos, 'bof');
end

function spkInd = spkIndex(post,ts)
M = length(ts);
spkInd = zeros(M,1);
for ii=1:M
    tdiff = (post-ts(ii)).^2;
    [m,ind] = min(tdiff);
    spkInd(ii) = ind;
end

%%%%%%%%%%%%% Head Direction Scripts %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decodes the target data.
function [dTargets,trackingColour] = decodeTargets(targets)

% Number of samples
numSamp = size(targets,2);

% Allocate memory to the array. 9 fields per sample: X-coord, Y-coord and
% 7 colour flag.
% Colour flag: 3=luminance, 4=rawRed, 5=rawGreen, 6=rawBlue, 7=pureRed,
% 8=pureGreen, 9=pureBlue.
dTargets = int16(zeros(numSamp,50,9));

for ii = 1:numSamp
    for jj = 1:50
        bitField = bitget(targets(jj,ii),1:32);
        if bitField(13)% Raw blue
            % Set the x-coord to the target
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            % Set the y-coord to the target
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,6) = 1;
        end
        if bitField(14) % Raw green
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,5) = 1;
        end
        if bitField(15) % Raw red
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,4) = 1;
        end
        if bitField(16) % Luminance
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,3) = 1;
        end
        if bitField(29) % Pure blue
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,9) = 1;
        end
        if bitField(30) % Puregreen
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,8) = 1;
        end
        if bitField(31) % Pure red
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,7) = 1;
        end
    end
end

% Find out what colours were used in the tracking
trackingColour = zeros(1,7);
if ~isempty(find(dTargets(:,:,3),1)) % Luminance
    trackingColour(1) = 1;
end
if ~isempty(find(dTargets(:,:,7),1)) % Pure Red
    trackingColour(2) = 1;
end
if ~isempty(find(dTargets(:,:,8),1)) % Pure Green
    trackingColour(3) = 1;
end
if ~isempty(find(dTargets(:,:,9),1)) % Pure Blue
    trackingColour(4) = 1;
end
if ~isempty(find(dTargets(:,:,4),1)) % Raw Red
    trackingColour(5) = 1;
end
if ~isempty(find(dTargets(:,:,5),1)) % Raw Green
    trackingColour(6) = 1;
end
if ~isempty(find(dTargets(:,:,6),1)) % Raw Blue
    trackingColour(7) = 1;
end

% Exctracts the individual coordinates for the centre of mass of each
% tracking diode. The red LEDs are assumed to be at the front and the green
% diodes are assumed to be at the back.
function [frontX,frontY,backX,backY] = extractPosition(targets,tracking)

ind = find(tracking(2:end));
if length(ind) <= 1
    % Need at least two colours to get head direction
    disp('ERROR: To few LED colours have been tracked. Not possible to find head direction')
    frontX = NaN;
    frontY = NaN;
    backX = NaN;
    backY = NaN;
    return
else
    if ~tracking(2) && ~tracking(5)
        disp('ERROR: Red LED has not been tracked')
        frontX = NaN;
        frontY = NaN;
        backX = NaN;
        backY = NaN;
        return
    end
    if ~tracking(3) && ~tracking(6)
        disp('ERROR: Green LED has not been tracked')
        frontX = NaN;
        frontY = NaN;
        backX = NaN;
        backY = NaN;
        return
    end
end

% Number of samples in the data
numSamp = size(targets,1);

% Allocate memory for the arrays
frontX = zeros(1,numSamp);
frontY = zeros(1,numSamp);
backX = zeros(1,numSamp);
backY = zeros(1,numSamp);

% Exctract the front coordinates (red LED)
if tracking(2) && ~tracking(5)
    % Pure red but not raw red
    for ii = 1:numSamp
        ind = find(targets(ii,:,7));
        if ~isempty(ind)
            frontX(ii) = mean(targets(ii,ind,1));
            frontY(ii) = mean(targets(ii,ind,2));
        end
    end
end
if ~tracking(2) && tracking(5)
    % Not pure red but raw red
    for ii = 1:numSamp
        ind = find(targets(ii,:,4));
        if ~isempty(ind)
            frontX(ii) = mean(targets(ii,ind,1));
            frontY(ii) = mean(targets(ii,ind,2));
        end
    end
end
if tracking(2) && tracking(5)
    % Both pure red and raw red
    for ii = 1:numSamp
        ind = find(targets(ii,:,7) | targets(ii,:,4));
        if ~isempty(ind)
            frontX(ii) = mean(targets(ii,ind,1));
            frontY(ii) = mean(targets(ii,ind,2));
        end
    end
end

% Exctract the back coordinates (green LED)
if tracking(3) && ~tracking(6)
    % Pure green but not raw green
    for ii = 1:numSamp
        ind = find(targets(ii,:,8));
        if ~isempty(ind)
            backX(ii) = mean(targets(ii,ind,1));
            backY(ii) = mean(targets(ii,ind,2));
        end
    end
end
if ~tracking(3) && tracking(6)
    % Not pure green but raw green
    for ii = 1:numSamp
        ind = find(targets(ii,:,5));
        if ~isempty(ind)
            backX(ii) = mean(targets(ii,ind,1));
            backY(ii) = mean(targets(ii,ind,2));
        end
    end
end
if tracking(3) && tracking(6)
    % Both pure green and raw green
    for ii = 1:numSamp
        ind = find(targets(ii,:,8) | targets(ii,:,5));
        if ~isempty(ind)
            backX(ii) = mean(targets(ii,ind,1));
            backY(ii) = mean(targets(ii,ind,2));
        end
    end
end

% Estimates lacking position samples using linear interpolation. When more
% than timeTreshold sec of data is missing in a row the data is left as
% missing.
%
% Raymond Skjerpeng 2006.
function [x] = interporPos(x,timeTreshold,sampRate)

% Turn off warning
warning('off','MATLAB:divideByZero');

% Number of samples that corresponds to the time threshold.
sampTreshold = floor(timeTreshold * sampRate);

% number of samples
numSamp = length(x);
% Find the indexes to the missing samples
temp1 = 1./x;
indt1 = isinf(temp1);
ind = indt1;
ind2 = find(ind==1);
% Number of missing samples
N = length(ind2);

if N == 0
    % No samples missing, and we return
    return
end

change = 0;

% Remove NaN in the start of the path
if ind2(1) == 1
    change = 1;
    count = 0;
    while 1
        count = count + 1;
        if ind(count)==0
            break
        end
    end
    x(1:count) = x(count);
end

% Remove NaN in the end of the path
if ind2(end) == numSamp
    change = 1;
    count = length(x);
    while 1
        count = count - 1;
        if ind(count)==0
            break
        end
    end
    x(count:numSamp) = x(count);
end

if change
    % Recalculate the missing samples
    temp1 = 1./x;
    indt1 = isinf(temp1);
    % Missing samples are where both x and y are equal to zero
    ind = indt1;
    ind2 = find(ind==1);
    % Number of samples missing
    N = length(ind2);
end

for ii = 1:N
    % Start of missing segment (may consist of only one sample)
    start = ind2(ii);
    % Find the number of samples missing in a row
    count = 0;
    while 1
        count = count + 1;
        if ind(start+count)==0
            break
        end
    end
    % Index to the next good sample
    stop = start+count;
    if start == stop
        % Only one sample missing. Setting it to the last known good
        % sample
        x(start) = x(start-1);
    else
        if count < sampTreshold
            % Last good position before lack of tracking
            x1 = x(start-1);
            % Next good position after lack of tracking
            x2 = x(stop);
            % Calculate the interpolated positions
            X = interp1([1,2],[x1,x2],1:1/count:2);
            % Switch the lacking positions with the estimated positions
            x(start:stop) = X;
            
            % Increment the counter (avoid estimating allready estimated
            % samples)
            ii = ii+count;
        else
            % To many samples missing in a row and they are left as missing
            ii = ii+count;
        end
    end
end


% Moving window mean smoothing filter
function [posx,posy] = posMeanFilter(posx,posy)

% Smooth samples with a mean filter over 15 samples
for cc = 8:length(posx)-7
    posx(cc) = nanmean(posx(cc-7:cc+7));
    posy(cc) = nanmean(posy(cc-7:cc+7));
end

% Calculates the direction of the head stage from the two set of
% coordinates. If one or both coordinate sets are missing for one samle the
% direction is set to NaN for that sample. Direction is also set to NaN for
% samples where the two coordinate set are identical. Returns the
% direction in degrees
function direct = headDirection(frontX,frontY,backX,backY)

% Number of position samples in data set
N = length(frontX);
direct = zeros(N,1);

for ii = 1:N
    
    if frontX(ii)==0 || backX(ii)==0
        % One or both coordinates are missing. No angle.
        direct(ii) = NaN;
        continue
    end
    
    % Calculate the difference between the coordinates
    xd = frontX(ii) - backX(ii);
    yd = frontY(ii) - backY(ii);
    
    if xd==0
        if yd==0
            % The two coordinates are at the same place and it is not
            % possible to calculate the angle
            direct(ii) = NaN;
            continue
        elseif yd>0
            direct(ii) = 90;
            continue
        else
            direct(ii) = 270;
            continue
        end
    end
    if yd==0
        if xd>0
            % Angle is zero
            continue
        else
            direct(ii) = 180;
            continue
        end
    end
    
    if frontX(ii)>backX(ii) && frontY(ii)>backY(ii)
        % Angle between 0 and 90 degrees
        direct(ii) = atan(yd/xd) * 360/(2*pi);
        
    elseif frontX(ii)<backX(ii) && frontY(ii)>backY(ii)
        % Angle between 90 and 180 degrees
        direct(ii) = 180 - atan(yd/abs(xd)) * 360/(2*pi);
        
    elseif frontX(ii)<backX(ii) && frontY(ii)<backY(ii)
        % Angle between 180 and 270 degrees
        direct(ii) = 180 + atan(abs(yd)/abs(xd)) * 360/(2*pi);
        
    else
        % Angle between 270 and 360 degrees
        direct(ii) = 360 - atan(abs(yd)/xd) * 360/(2*pi);
    end
end