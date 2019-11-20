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

function MCWPLI_opto(inFile,numperms)  % width - number of cycles in wavelet (> 5 advisable); order - whitening

img_text = 'on';

fid = fopen(inFile,'r');
if fid == -1
    msgbox('Could not open the input file! Make sure the filname and path are correct.','ERROR');
end

% Get sessions and csc-file list from input file
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

%make directory to save figures
savedir = cd;
dirInfo = dir(cd);
found = 0;
for kk=1:size(dirInfo,1)
    if dirInfo(kk).isdir
        if strcmp(dirInfo(kk).name,strcat('MCWPLI_opto','\'))
            found = 1;
        end
    end
end
if found==0
    mkdir(savedir,strcat('MCWPLI_opto','\'));
end

veloff = [];
velon = [];
ImXoff = [];
ImXon = [];

for jj=1:numsessions
    
    % Get position data
    % Set the field selection for reading the video files. 1 = Add parameter, 0 = skip
    % parameter

        
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
        file = strcat(sessions{jj},'vt1.nvt');
        [t{jj}, x, y] = Nlx2MatVT(file,fieldSelection,extractHeader,extractMode);
        ind = find(x == 0);
        
        t{jj}(ind) = [];
        x(ind) = [];
        y(ind) = [];
        tidxoff{jj} = [];
        tidxon{jj} = [];
        
        % Adjust input data
        D = 95;
        xscale = D/(max(x)-min(x));
        yscale = D/(max(y)-min(y));
        x = x * xscale;
        y = y * yscale;
        
        % Convert timestamps to seconds
        t{jj} = t{jj}/1000000;
        % Do median filtering to supress off-values
        %[x, y, t{jj}] = medianFilter(x, y, t{jj});
        % Smoothing the position samples with a simple mean filter
        for cc=8:length(x)-7
            x(cc) = nanmean(x(cc-7:cc+7));
            y(cc) = nanmean(y(cc-7:cc+7));
        end
        
        % With many off-values in a row the median filter doesn't remove all, must
        % apply an additional off-value filter.
        %[x, y, t{jj}] = offValueFilter(x, y, t{jj});
        
        [x,y] = centreBox(x,y);
        x_cen = (max(x) - min(x))/2 + min(x);
        y_cen = (max(y) - min(y))/2 + min(y);
        x = x - x_cen;
        y = y - y_cen;
        
        %compute velocity
        vel{jj} = zeros(1,length(x));
        for v = 2:1:(length(x)-1)
            vel{jj}(v) = sqrt((x(v+1)-x(v-1))^2 + (y(v+1)-y(v-1))^2)/(t{jj}(v+1)-t{jj}(v-1));
        end
        vel{jj}(end) = sqrt((x(end)-x(end-1))^2 + (y(end)-y(end-1))^2)/(t{jj}(length(x))-t{jj}(length(x)-1));
        vel{jj}(vel{jj}>=40) = 0.5*(vel{jj}(circshift((vel{jj}>=40),[-3,0]))+vel{jj}(circshift((vel{jj}>=40),[3,0])));
        vel{jj} = smooth(vel{jj},15);
        
        %Set all vars to common time
        mt{jj} = min(t{jj});
        t{jj} = t{jj} - mt{jj};
        
        if jj == 1 || jj == 2
            efile = strcat(sessions{jj},'Events.nev');
            FS = [1 0 1 0 1]; EH = 0; EM = 1;
            [TimeStamps, TTLs, EventStrings] = Nlx2MatEV(efile,FS,EH,EM);
            TTLs(1) = []; TimeStamps(1) = []; TimeStamps = TimeStamps/1000000;
            if isempty(TTLs)
                msgbox('No TTLs in input file!','ERROR');
            end
            TS = TimeStamps - mt{jj};
            
            
            for ts=2:length(TS)
                tidx = find(t{jj}>TS(ts-1)&t{jj}<=TS(ts));
                if TTLs(ts) == 0
                    tidxoff{jj} = [tidxoff{jj} tidx];
                else
                    tidxon{jj} = [tidxon{jj} tidx];
                end
            end
            velon = [velon vel{jj}(tidxon{jj})'];
            veloff = [veloff vel{jj}(tidxoff{jj})'];
        end
    
    if jj == 1 || jj == 2        
        
        %resize ImX to one sample per vel
        data = load([sessions{jj},strcat('MCWPLI','\MCWPLI'),'.mat']);
        tt = data.tt - mt{jj};
        ImX_tt = data.MCImX;
        freqVec = data.freq;
        clear data
        ImX_t = zeros(length(freqVec),length(t));
        w = median(diff(t{jj}));
        for i = 1:length(t{jj})
            ImX_t(:,i) = nanmean(ImX_tt(:,abs(tt-t{jj}(i))<=w),2);
        end
        clear ImX_tt
        
        ImXon = [ImXon ImX_t(:,tidxon{jj})];
        ImXoff = [ImXoff ImX_t(:,tidxoff{jj})];
        
        clear ImX_t
    end
end

wpli(:,1) = abs(nansum(ImXoff,2))./nansum(abs(ImXoff),2);
wpli(:,2) = abs(nansum(ImXon,2))./nansum(abs(ImXon),2);

for jj=3:numsessions
    
    %resize ImX to one sample per vel
    data = load([sessions{jj},strcat('MCWPLI','\MCWPLI'),'.mat']);
    tt = data.tt - mt{jj};
    ImX_tt = data.MCImX;
    freqVec = data.freq;
    clear data
    ImX_t{jj} = zeros(length(freqVec),length(t));
    w = median(diff(t{jj}));
    for i = 1:length(t{jj})
        ImX_t{jj}(:,i) = nanmean(ImX_tt(:,abs(tt-t{jj}(i))<=w),2);
    end
    clear ImX_tt
    
end

firstlight = linspace(0,30,numperms);
p_onwpli = zeros(length(freqVec),numperms);
p_offwpli = zeros(length(freqVec),numperms);

for b = 1:numperms
    ImXon = []; ImXoff = [];
    velon = []; veloff = [];
    
    for jj = 3:numsessions
        tidxon{jj} = []; tidxoff{jj} = [];
        TS = firstlight(b):30:max(t{jj});
        pTTLs = ones(size(TS));
        pTTLs(1:2:end) = 0;
        
        for ts=2:length(TS)
            tidx = find(t{jj}>TS(ts-1)&t{jj}<=TS(ts));
            if pTTLs(ts) == 0
                tidxoff{jj} = [tidxoff{jj} tidx];
            else
                tidxon{jj} = [tidxon{jj} tidx];
            end
        end
        
        velon = [velon vel{jj}(tidxon{jj})'];
        veloff = [veloff vel{jj}(tidxoff{jj})'];
        ImXon = [ImXon ImX_t{jj}(:,tidxon{jj})];
        ImXoff = [ImXoff ImX_t{jj}(:,tidxoff{jj})];
        
    end
    
    p_offwpli(:,b) = abs(nansum(ImXoff,2))./nansum(abs(ImXoff),2);
    p_onwpli(:,b) = abs(nansum(ImXon,2))./nansum(abs(ImXon),2);
end

pstd = std([p_onwpli-p_offwpli p_offwpli-p_onwpli],0,2);

close all
figure(1);
subplot(1,2,1)
[ci] = prctile([(p_onwpli-p_offwpli)./p_offwpli*100 (p_offwpli-p_onwpli)./p_onwpli*100]',[2.5 97.5]);
shadedplot(freqVec, ci(2,:), ci(1,:),[1 0 0],[1 0.25 0.25]); alpha(0.5);hold on
plot(freqVec,mean([p_onwpli-p_offwpli p_offwpli-p_onwpli],2),'-r')
plot(freqVec,(wpli(:,2)-wpli(:,1))./wpli(:,1)*100,'-k')
xlabel('frequency (Hz)')
ylabel('percent change from no light')

subplot(1,2,2)
[ci] = prctile([p_onwpli-p_offwpli p_offwpli-p_onwpli]'./repmat(pstd,1,2*numperms)',[2.5 97.5]);
shadedplot(freqVec, ci(2,:), ci(1,:),[1 0 0],[1 0.25 0.25]); alpha(0.5);hold on
plot(freqVec,mean([p_onwpli-p_offwpli p_offwpli-p_onwpli]'./repmat(pstd,1,2*numperms)'),'-r')
plot(freqVec,(wpli(:,2)-wpli(:,1))./pstd,'-k')
xlabel('frequency (Hz)')
ylabel('z-score from control')

figImage = sprintf('%s%s%s%s',savedir,strcat('\MCWPLI_opto','\'),'MCWPLI_opto.fig');
bmpImage = sprintf('%s%s%s%s',savedir,strcat('\MCWPLI_opto','\'),'MCWPLI_opto.bmp');
saveas(gcf,figImage,'fig');
saveas(gcf,bmpImage,'bmp');



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




