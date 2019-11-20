function [LFPmaps] = LFPmap_batch(inFile,freqVec,startmin,endmin,width) 

h = 4; % Smoothing factor when calculating the LFPmap
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

% Read the input data
timeMaps = cell(numsessions,1);
P = cell(numsessions,numchannels);
LFPmaps = cell(numsessions,numchannels,length(freqVec));
% fam_dis_idx = zeros(numsessions,numchannels,length(freqVec));
% test_dis_idx = zeros(numsessions,numchannels,length(freqVec));
% pathCoord = cell(numsessions,3);

for ii = 1:numsessions
    disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));
    % Check if subdir for storing images are present. If not, it is
    % created
    dirInfo = dir(sessions{ii});
    found = 0;
    foundv = 0;
    foundt = 0;
    for kk=1:size(dirInfo,1)
        if dirInfo(kk).isdir
            if strcmp(dirInfo(kk).name,strcat('LFPmaps','\'))
                found = 1;
            end
            if strcmp(dirInfo(kk).name,strcat('Vmaps','\'))
                foundv = 1;
            end
            if strcmp(dirInfo(kk).name,strcat('Tmaps','\'))
                foundt = 1;
            end
        end
    end
    if found==0
        mkdir(sessions{ii},strcat('LFPmaps','\'));
    end
    if foundv==0
        mkdir(sessions{ii},strcat('Vmaps','\'));
    end
    if foundt==0
        mkdir(sessions{ii},strcat('Tmaps','\'));
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
    xmin = xmin-xcen;
    xmax = xmax-xcen;
    ymin = ymin-ycen;
    ymax = ymax-ycen;
    
    % Convert timestamps to seconds
    t = t/1000000;
    
    sLength = xmax-xmin;
    bins = 30;
    binWidth = sLength/bins;
    mapAxis = (-sLength/2+binWidth/2):binWidth:(sLength/2-binWidth/2);     
   
    x(t-min(t)>60*(endmin)) = [];
    y(t-min(t)>60*(endmin)) = [];
    t(t-min(t)>60*(endmin)) = [];
    x(t-min(t)<60*(startmin)) = [];
    y(t-min(t)<60*(startmin)) = [];
    t(t-min(t)<60*(startmin)) = [];
    
    %compute velocity and velocity map
    vel = zeros(1,length(x));
    for v = 2:1:(length(x)-1)
        vel(v) = sqrt((x(v+1)-x(v-1))^2 + (y(v+1)-y(v-1))^2)/(t(v+1)-t(v-1));
    end
    vel(end) = sqrt((x(end)-x(end-1))^2 + (y(end)-y(end-1))^2)/(t(length(x))-t(length(x)-1));
    vel(vel>=40) = 0.5*(vel(circshift((vel>=40),-3))+vel(circshift((vel>=40),3)));
    vel = smooth(vel,15);  

    % Calculate the time map for this session
    timeMaps{ii,1} = findTimeMap(x,y,t,mapAxis);
    dt = diff(t); dt(end+1) = dt(end);
    Tmap = LFPmap(dt,x',y',h,mapAxis);
    Tmap(timeMaps{ii}==0) = NaN;
    Vmap = LFPmap(vel',x',y',h,mapAxis);    
    Vmap(timeMaps{ii}==0) = NaN;
    
    %plot and save velocity map
    colormap(hot(500));
    imagesc(mapAxis,mapAxis,Vmap,[min(min(Vmap)) max(max(Vmap))]); axis xy; colorbar;
    colordata = colormap; colordata(1,:) = [0.6 0.8 1];colormap(colordata);
    c = colorbar; ylabel(c,'running speed (cm/sec)');
    xlim([xmin xmax]); ylim([ymin ymax]);
    vlineX = [(xmin + (xmax-xmin)/2) (xmin + (xmax-xmin)/2)];
    vlineY = [ymax ymin];
    hlineX = [xmin xmax];
    hlineY = [(ymin + (ymax-ymin)/2) (ymin + (ymax-ymin)/2)];
    hold on; plot(x,y,'Color',[0 0 0]);
    plot (vlineX,vlineY, 'k');
    plot (hlineX,hlineY, 'k'); hold off;
    title(['Running Speed Map: ',num2str(startmin),'-',num2str(endmin),]);
    figImage = sprintf('%s%s%s%s',sessions{ii},strcat('Vmaps','\Vmap_'),num2str(startmin),'-',num2str(endmin),'10.fig');
    bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('Vmaps','\Vmap_'),num2str(startmin),'-',num2str(endmin),'10.bmp');
    %epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('Vmaps','\Vmap_'),num2str(startmin),'-',num2str(endmin),'.eps');
    saveas(gcf,figImage,'fig');
    saveas(gcf,bmpImage,'bmp');
    %saveas(gcf,epsImage,'eps');  
    close all
    colormap(hot(500));
    imagesc(mapAxis,mapAxis,Tmap,[min(min(Tmap)) max(max(Tmap))]); axis xy; colorbar;
    c = colorbar; ylabel(c,'time (sec)')
    colordata = colormap; colordata(1,:) = [0.6 0.8 1];colormap(colordata);
    xlim([xmin xmax]); ylim([ymin ymax]);    
    vlineX = [(xmin + (xmax-xmin)/2) (xmin + (xmax-xmin)/2)];
    vlineY = [ymax ymin];
    hlineX = [xmin xmax];
    hlineY = [(ymin + (ymax-ymin)/2) (ymin + (ymax-ymin)/2)];
    hold on; plot(x,y,'Color',[0 0 0]);
    plot (vlineX,vlineY, 'k');
    plot (hlineX,hlineY, 'k'); hold off;
    title(['Time Map: ',num2str(startmin),'-',num2str(endmin),]);
    figImage = sprintf('%s%s%s%s',sessions{ii},strcat('Tmaps','\Tmap_'),num2str(startmin),'-',num2str(endmin),'10.fig');
    bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('Tmaps','\Tmap_'),num2str(startmin),'-',num2str(endmin),'10.bmp');
    %epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('Tmaps','\Tmap_'),num2str(startmin),'-',num2str(endmin),'.eps');
    saveas(gcf,figImage,'fig');
    saveas(gcf,bmpImage,'bmp');
    %saveas(gcf,epsImage,'eps');  
    
    % END COMPUTING POSITION DATA
    
    % Load data from the .ncs files, make plots, and store them
    for jj=1:numchannels
        disp('Make plots and store them to files');
        disp(sprintf('%s%i',' CSC ',jj, ' of ',numchannels));
        file = [sessions{ii},channels{jj}];
        [samples,ts,tt, Fs, bv, ir] = loadEEG2(file);
        ch_X = bv*samples; ts = tt;
        ch_X(ts-min(ts)>60*endmin) = [];
        ch_X(ts-min(ts)<60*startmin) = [];   
        ts(ts-min(ts)>60*endmin) = []; 
        ts(ts-min(ts)<60*startmin) = [];               
        [P{ii,jj}] = traces2TFR([ch_X ch_X],freqVec,Fs,width);      
        P{ii,jj} = mean(P{ii,jj});
        P{ii,jj} = zscore(P{ii,jj});
        
        %Get position of P 
        [Px,Py] = PPos(ts,x,y,t);
        %Calculate LFP map
        LFPmaps{ii,jj} = LFPmap(P{ii,jj},Px,Py,h,mapAxis);
        % Remove unvisited parts of the box from the map
        if length(LFPmaps{ii,jj})>1
            LFPmaps{ii,jj}(timeMaps{ii}==0) = NaN;
        end
        
        data.P = P{ii,jj}; data.freq = freqVec; data.LFPmaps = LFPmaps(ii,jj,:);   %save data as struct
        data.Px = Px; data.Py = Py; data.ts = ts; data.x = x; data.y = y; data.t = t;
        filename = sprintf('%s%s%s%s',sessions{ii},strcat('LFPmaps','\'),channels{jj}(1:end-4),'_',num2str(min(freqVec)),'-',num2str(max(freqVec)),' Hz_cont','.mat');
        save(filename,'-struct','data');        
        disp('Plot maps and store them to files');    
        figure(1);colormap(hot(500));
        imagesc(mapAxis,mapAxis,LFPmaps{ii,jj},[min(min(LFPmaps{ii,jj})) max(max(LFPmaps{ii,jj}))]); axis xy; colorbar;
        c = colorbar; ylabel(c,'z-scored power');
        colordata = colormap; colordata(1,:) = [0.6 0.8 1];colormap(colordata);
        xlim([xmin xmax]); ylim([ymin ymax]);
        vlineX = [(xmin + (xmax-xmin)/2) (xmin + (xmax-xmin)/2)];
        vlineY = [ymax ymin];
        hlineX = [xmin xmax];
        hlineY = [(ymin + (ymax-ymin)/2) (ymin + (ymax-ymin)/2)];
        hold on; plot(x,y,'Color',[0 0 0]); 
        plot (vlineX,vlineY, 'k');
        plot (hlineX,hlineY, 'k'); hold off;
        title(['LFP Map: ',num2str(min(freqVec)),'-',num2str(max(freqVec)),' Hz: ',file]);
        figImage = sprintf('%s%s%s%s',sessions{ii},strcat('LFPmaps','\LFPmap_'),channels{jj}(1:end-4),'_',num2str(min(freqVec)),'-',num2str(max(freqVec)),' Hz','10.fig');
        bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('LFPmaps','\LFPmap_'),channels{jj}(1:end-4),'_',num2str(min(freqVec)),'-',num2str(max(freqVec)),' Hz','10.bmp');
%         epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('LFPmaps','\LFPmap_'),channels{jj}(1:end-4),'_',num2str(min(freqVec)),'-',num2str(max(freqVec)),' Hz','.10eps');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');
%         saveas(gcf,epsImage,'eps');

        for f = 1:length(freqVec)
%             LFPmaps{ii,jj,f} = reshape(tiedrank(reshape(LFPmaps{ii,jj,f},1,60*60)),60,60);
%             novel = reshape(LFPmaps{ii,jj,f}(:,31:60),1,60*30);
%             %novel = reshape(LFPmaps{ii,jj}(:,1:30),1,60*30);
%             novel(isnan(novel)) = [];
%             familiar = reshape(LFPmaps{ii,jj,f}(:,1:30),1,60*30);
%             %familiar = reshape(LFPmaps{ii,jj}(:,31:60),1,60*30);
%             familiar(isnan(familiar)) = [];
%             fam_dis_idx(ii,jj,f) = median(novel)/(median(novel)+median(familiar));
            
%             novel = reshape(LFPmaps{ii,jj,f}(1:30,1:30),1,30*30);
%             %novel = reshape(LFPmaps{ii,jj}(1:30,31:60),1,30*30);
%             novel(isnan(novel)) = [];
%             familiar = reshape(LFPmaps{ii,jj,f}(31:60,31:60),1,30*30);
%             %familiar = reshape(LFPmaps{ii,jj}(1:30,1:30),1,30*30);
%             familiar(isnan(familiar)) = [];
%             test_dis_idx(ii,jj,f) = median(novel)/(median(novel)+median(familiar));
        end

    end
end

% Discrim_idx.fam = fam_dis_idx;
% Discrim_idx.test = test_dis_idx;
% 
% filename = sprintf('%s%s%s%s','Discrim_idx_',num2str(min(freqVec)),'-',num2str(max(freqVec)),' Hz_cont','.mat');
% save(filename,'-struct','Discrim_idx');


%__________________________________________________________________________
%
% Field functions
%__________________________________________________________________________

% timeMap will find the amount of time the rat spends in each bin of the
% box. The box is divided into the same number of bins as in the rate map.
% posx, posy, post: Path coordinates and time stamps.
% bins: Number of bins the box is diveded int (bins x bins)
% extremal: minimum and maximum positions.
function timeMap = findTimeMap(posx,posy,post,mapAxis)

% Duration of trial
duration = post(end)-post(1);
% Average duration of each position sample
sampDur = duration/length(posx);
timeMap = zeros(length(mapAxis),length(mapAxis));
% Find number of position samples in each bin
bwidth = mode(diff(mapAxis));
yy = 0;
for ii = mapAxis
    yy = yy + 1;
    xx = 0;
    aty = posy >= ii-bwidth/2 & posy <= ii+bwidth/2;
    for jj = mapAxis
        xx = xx + 1;
        atx = posx >= jj-bwidth/2 & posx <= jj+bwidth/2;
        timeMap(yy,xx) = sum(atx & aty);
    end
end
% Convert to time spent in bin
timeMap = timeMap * sampDur;

% Calculates the rate map.
function map = LFPmap(P,Px,Py,h,mapAxis)
invh = 1/h;
map = zeros(length(mapAxis),length(mapAxis));
yy = 0;
for j = mapAxis
    yy = yy + 1;
    xx = 0;
    for i = mapAxis
        xx = xx + 1;
        map(yy,xx) = P_estimator(P,Px,Py,i,j,invh);
    end
end

% Estimate the mean P for one position value
function p = P_estimator(P,Px,Py,i,j,invh)
% edge-corrected kernel density estimator
conv_sum = P*gaussian_kernel((Px-i),(Py-j),invh);
edge_corrector =  sum(gaussian_kernel((Px-i),(Py-j),invh));
p = conv_sum/edge_corrector; % regularised firing rate for "wellbehavedness"                                                      
if P(end) == P(end-1) && mean(P)>0.02 && mean(P) < 0.05 %i.e. if time map
    p=conv_sum;
end
                                                    
% Gaussian kernel for the rate calculation
function r = gaussian_kernel(x,y,invh)
r =  invh^2/sqrt(2*pi)*exp(-0.5*(x.*x*invh^2 + y.*y*invh^2));

% Finds the position to the spikes
function [Px,Py] = PPos(ts,x,y,t)
N = length(ts);
Px = zeros(N,1);
Py = zeros(N,1);
for ii = 1:N
    tdiff = (t-ts(ii)).^2;
    [m,ind] = min(tdiff);
    Px(ii) = x(ind(1));
    Py(ii) = y(ind(1));
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

function drawfield(map,faxis,cmap,maxrate)
   
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
%    if strcmp(text,'on')
%        title(strcat(cellid,'¤(0 - ',num2str(maxrate), ' Hz)','¤',cell_file),'FontSize',20);
%    end


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

