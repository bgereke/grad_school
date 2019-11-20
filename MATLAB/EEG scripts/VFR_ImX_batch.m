% CFC('inFile.txt')
%
% CFC generates cross frequency coherences, plots them and stores them to images and .mat files.
% Multiple sessions will be read from the CSC's specififed in
% 'CSCList.txt'. 'TTList.txt' is necessary for spike data scripts that use
% the same 'infile.txt'.
%
% The input file must be on the following format.
%
% C:\Data\TTList.txt
% C:\Data\CSCList.txt
% C:\Data\Begin 1
% C:\Data\Begin 2
% C:\Data\Begin 3
% C:\Data\Begin 4
% and so on ...
%
% 'CSCList.txt' contains a list of the Neuralynx .csc files to be analyzed.
% All plots will be stored to both bmp and eps imagefiles to a subdirectory in
% the data folder called CFC_plots.

function [VFR] = VFR_ImX_batch(inFile,freqVec)

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
    if ii == 0
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

% read the file names from the csc-file list
cscid = fopen(cscList,'r');
jj = 1;
while ~feof(cscid)
       str = fgetl(cscid);
       channels(jj) = {str};
       jj = jj+1;
end
numchannels = jj-1;
cscid = fclose('all');

% read the file names of the references for the channels to be used for
% averaging
refid = fopen(refList,'r');
jj = 1;
while ~feof(refid)
       str = fgetl(refid);
       refs(jj) = {str};
       jj = jj+1;
end

% read the file names of the channels to be used for averaging
avgid = fopen(ch4avg,'r');
jj = 1;
while ~feof(avgid)
       str = fgetl(avgid);
       avgs(jj) = {str};
       jj = jj+1;
end
numrefs = jj-1;

for ii = 1:numsessions
    disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));
    % Check if subdir for storing images are present. If not, it is
    % created
    dirInfo = dir(sessions{ii});
    found = 0;
    for kk=1:size(dirInfo,1)
        if dirInfo(kk).isdir
            if strcmp(dirInfo(kk).name,strcat('VFR_ImX','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(sessions{ii},strcat('VFR_ImX','\'));
    end
    
    % Get position data
    fieldSelection(1) = 1; % Timestamps
    fieldSelection(2) = 1; % Extracted X
    fieldSelection(3) = 1; % Extracted Y
    fieldSelection(4) = 0; % Extracted Angel
    fieldSelection(5) = 0; % Targets
    fieldSelection(6) = 0; % Points
    % Do we return header 1 = Yes, 0 = No.
    extractHeader = 0;
    % 5 different extraction modes, see help file for Nlx2MatVt
    extractMode = 1; % Extract all data
    file = strcat(sessions{ii},'vt1.nvt');
    [t, x, y] = Nlx2MatVT(file,fieldSelection,extractHeader,extractMode);  
    
    ind = find(x == 0);
    t(ind) = [];
    x(ind) = [];
    y(ind) = [];
    % Do median filtering to supress off-values
    %[x, y, t] = medianFilter(x, y, t);
    % Smoothing the position samples with a simple mean filter
    for cc=8:length(x)-7
        x(cc) = nanmean(x(cc-7:cc+7));
        y(cc) = nanmean(y(cc-7:cc+7));
    end
    y = -y + max(y)-min(y) + 2*max(y); % reflects posisitons so they are consistent with orientation in recording room
    [hx,xc] = hist(x,50); [hy,yc] = hist(y,50);
    th = 50;
    xmin = xc(find(hx>th,1,'first'));
    xmax = xc(find(hx>th,1,'last'));
    ymin = yc(find(hy>th,1,'first'));
    ymax = yc(find(hy>th,1,'last'));

    % Adjust input data
    scalex = 40/(xmax-xmin);
    scaley = 40/(ymax-ymin);
    x = x * scalex;
    y = y * scaley;
    xmin = xmin * scalex;
    ymin = ymin * scaley;
    xmax = xmax * scalex;
    ymax = ymax * scaley;
    xcen = xmin+(xmax-xmin)/2;
    ycen = ymin+(ymax-ymin)/2;
    x = x - xcen;
    y = y - ycen;
    
    % Convert timestamps to seconds
    t = t/1000000;
    
    %compute velocity and velocity map
    vel = zeros(1,length(x));
    for v = 2:1:(length(x)-1)
        vel(v) = sqrt((x(v+1)-x(v-1))^2 + (y(v+1)-y(v-1))^2)/(t(v+1)-t(v-1));
    end
    vel(end) = sqrt((x(end)-x(end-1))^2 + (y(end)-y(end-1))^2)/(t(length(x))-t(length(x)-1));
    vel(vel>=40) = 0.5*(vel(circshift((vel>=40),-3))+vel(circshift((vel>=40),3)));
    vel = smooth(vel,15);
    
%      % Make Average Reference
%     avref = 0;
%     for jj=1:numrefs
%         file = [sessions{ii},avgs{jj}];
%         [samples,ts,tt, Fs, bv, ir] = loadEEG2(file);
%         ch_X = bv*samples;
%         cscnum = avgs{jj}(4:end-4); % get  tt# corresponding to avgsjj 
%         cscnum = str2num(cscnum);   
%         ttnum = ceil(cscnum/4);
%         if ~strcmp(refs{ttnum},'G')
%             rfile = [sessions{ii},refs{ttnum}];
%             [samples,ts,tt, Fs, bv, ir] = loadEEG2(rfile);
%             ch_X = ch_X + bv*samples;
%         end
%         avref = avref + ch_X/numrefs;
%         clear samples ch_X
%     end
    
    % Load data from the .ncs files, make plots, and store them
    
    disp('Make plots and store them to files'); 
    xfile = [sessions{ii},channels{1}];
    [samples,ts,tt, Fs, bv, ir] = loadEEG2(xfile);
    x = bv*samples;
    yfile = [sessions{ii},channels{2}];
    [samples,ts,tt, Fs, bv, ir] = loadEEG2(yfile);
    y = bv*samples;
    %set recording channel against ground and then against average reference
    xcscnum = channels{1}(4:end-4); % get  tt# corresponding to avgsjj
    xcscnum = str2num(xcscnum);
    xttnum = ceil(xcscnum/4);
    ycscnum = channels{1}(4:end-4); % get  tt# corresponding to avgsjj
    ycscnum = str2num(ycscnum);
    yttnum = ceil(ycscnum/4);
    if ~strcmp(refs{xttnum},'G')
        xrfile = [sessions{ii},refs{xttnum}];
        [samples,ts,tt, Fs, bv, ir] = loadEEG2(xrfile);
        x = x + bv*samples;
    end
    if ~strcmp(refs{yttnum},'G')
        yrfile = [sessions{ii},refs{yttnum}];
        [samples,ts,tt, Fs, bv, ir] = loadEEG2(yrfile);
        y = y + bv*samples;
    end
%     x = x - avref;
%     y = y - avref;

    %resize vel to one sample per tt
    v = zeros(1,length(tt));
    [ttmin, ttidx] = min((tt-t(1)).^2);
    v(1:ttidx) = vel(1);
    old = ttidx+1;
    for i = 2:length(t)
        [ttmin, ttidx] = min((tt-t(i)).^2);
        v(old:ttidx) = linspace(vel(i-1),vel(i),ttidx-old+1);
        old = ttidx+1;
    end
    if old < length(tt)
        v(old:end) = vel(end);
    end
    vel = v;
    vbins = min(vel)+0.1:0.8*max(vel);
    
    [VFR] = VFR_ImX(x,y,freqVec,6,Fs,vel,vbins);
    VFR = VFR;
    data.freq = freqVec; data.vbins = vbins; data.vfr = VFR;   %save data as struct
    filename = sprintf('%s%s%s%s',sessions{ii},strcat('VFR_ImX','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'.mat');
    save(filename,'-struct','data');
    figure(1);
    uimagesc(vbins(1:end-1),freqVec,VFR); axis xy; colorbar
    xlabel('log-running speed');ylabel('Frequency (Hz)');
    set(gca,'xtick',[5 10 20])
    set(gca,'XTickLabel',[5 10 20])
    set(gca,'ytick',[20 40 80])
    set(gca,'YTickLabel',[20 40 80])
%     imagesc(vbins(1:end-1),freqVec,VFR); axis xy; colorbar
    xlabel('running speed (cm/sec)');ylabel('frequency (Hz)');
    title(strcat('Average Running Speed Triggered ImX: ',channels{1}(1:end-4),'x',channels{2}(1:end-4)));
    figImage = sprintf('%s%s%s%s',sessions{ii},strcat('VFR_ImX','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'.fig');
    bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('VFR_ImX','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'.bmp');
    epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('VFR_ImX','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'.eps');
    saveas(gcf,figImage,'fig');
    saveas(gcf,bmpImage,'bmp');
    saveas(gcf,epsImage,'eps');
    
end

%_______________________________________________________________________
%other functions
%_______________________________________________________________________

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



