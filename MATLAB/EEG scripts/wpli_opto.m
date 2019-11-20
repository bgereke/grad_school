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

function wpli_opto(inFile,freqVec,width,numperms)  % width - number of cycles in wavelet (> 5 advisable); order - whitening

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

%make directory to save figures
savedir = cd;
dirInfo = dir(cd);
found = 0;
for kk=1:size(dirInfo,1)
    if dirInfo(kk).isdir
        if strcmp(dirInfo(kk).name,strcat('OptoWPLI','\'))
            found = 1;
        end
    end
end
if found==0
    mkdir(savedir,strcat('OptoWPLI','\'));
end

veloff = [];
velon = [];

for ii = 1:numchannels-1
    for kk = ii+1:numchannels
        
        
        %declare main vars        
        ImXoff = [];
        ImXon = [];
        ImXoffss = [];
        ImXonss = [];
        wpli  = zeros(length(freqVec),2);
        outsum   = zeros(length(freqVec),2);
        outsumW  = zeros(length(freqVec),2);
        
        disp('Make plots and store them to files');
        %disp(sprintf('%s%i',' CSC ',ii, ' of ',numchannels));
        
        for jj=1:numsessions-2
            
            if kk == 2
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
                [t, x, y] = Nlx2MatVT(file,fieldSelection,extractHeader,extractMode);
                ind = find(x == 0);
                
                t(ind) = [];
                x(ind) = [];
                y(ind) = [];
                
                % Adjust input data
                D = 95;
                xscale = D/(max(x)-min(x));
                yscale = D/(max(y)-min(y));
                x = x * xscale;
                y = y * yscale;
                
                % Convert timestamps to seconds
                t = t/1000000;
                % Do median filtering to supress off-values
                %[x, y, t] = medianFilter(x, y, t);
                % Smoothing the position samples with a simple mean filter
                for cc=8:length(x)-7
                    x(cc) = nanmean(x(cc-7:cc+7));
                    y(cc) = nanmean(y(cc-7:cc+7));
                end
                
                % With many off-values in a row the median filter doesn't remove all, must
                % apply an additional off-value filter.
                %[x, y, t] = offValueFilter(x, y, t);
                
                [x,y] = centreBox(x,y);
                x_cen = (max(x) - min(x))/2 + min(x);
                y_cen = (max(y) - min(y))/2 + min(y);
                x = x - x_cen;
                y = y - y_cen;
                
                %compute velocity
                vel = zeros(1,length(x));
                for v = 2:1:(length(x)-1)
                    vel(v) = sqrt((x(v+1)-x(v-1))^2 + (y(v+1)-y(v-1))^2)/(t(v+1)-t(v-1));
                end
                vel(end) = sqrt((x(end)-x(end-1))^2 + (y(end)-y(end-1))^2)/(t(length(x))-t(length(x)-1));
                vel(vel>=40) = 0.5*(vel(circshift((vel>=40),-3))+vel(circshift((vel>=40),3)));
                vel = smooth(vel,15);
            end
            
            efile = strcat(sessions{jj},'Events.nev');
            FS = [1 0 1 0 1]; EH = 0; EM = 1;
            [TimeStamps, TTLs, EventStrings] = Nlx2MatEV(efile,FS,EH,EM);
            TTLs(1) = []; TimeStamps(1) = []; TimeStamps = TimeStamps/1000000;
            if length(TTLs) == 0
                msgbox('No TTLs in input file!','ERROR');
            end
            
            
            filex = [sessions{jj},channels{ii}];
            [samples,ts,tt, Fs, bv, ir] = loadEEG2(filex);
            x = bv*samples;
            filey = [sessions{jj},channels{kk}];
            [samples,ts,tt, Fs, bv, ir] = loadEEG2(filey);
            y = bv*samples;
            %set recording channel against ground and then against average reference
            xcscnum = channels{1}(4:end-4); % get  tt# corresponding to avgsjj
            xcscnum = str2num(xcscnum);
            xttnum = ceil(xcscnum/4);
            ycscnum = channels{2}(4:end-4); % get  tt# corresponding to avgsjj
            ycscnum = str2num(ycscnum);
            yttnum = ceil(ycscnum/4);
            %         if ~strcmp(refs{xttnum},'G')
            %             xrfile = [sessions{jj},refs{xttnum}];
            %             [samples,ts,tt, Fs, bv, ir] = loadEEG2(xrfile);
            %             x = x + bv*samples;
            %         end
            %         if ~strcmp(refs{yttnum},'G')
            %             yrfile = [sessions{jj},refs{yttnum}];
            %             [samples,ts,tt, Fs, bv, ir] = loadEEG2(yrfile);
            %             y = y + bv*samples;
            %         end
            %     x = x - avref;
            %     y = y - avref;
            clear samples
            
            [ImX] = traces2ImX(x,y,freqVec,Fs,width);
            
            %Set all vars to common time
            ImX(:,tt<TimeStamps(1))=[];tt(tt<TimeStamps(1))=[];
            vel(t<TimeStamps(1)) = []; t(t<TimeStamps(1))=[];
            ImX(:,tt>TimeStamps(end))=[];tt(tt>TimeStamps(end))=[];
            vel(t>TimeStamps(end)) = [];t(t>TimeStamps(end))=[];
            tt = tt - min(TimeStamps);
            t = t - min(TimeStamps);
            TS = TimeStamps - min(TimeStamps);
            
            if kk == 2
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
            end
            
            for ts=2:length(TS)
                tidx = tt>TS(ts-1)&tt<=TS(ts);
                if TTLs(ts) == 0
                    if kk == 2
                        veloff = [veloff vel(tidx)];
                    end
                    ImXoff = [ImXoff ImX(:,tidx)];
                else
                    if kk == 2
                        velon = [velon vel(tidx)];
                    end
                    ImXon = [ImXon ImX(:,tidx)];
                end
            end
        end
        for jj=1:numsessions-2
            
            if kk == 2
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
                [t, x, y] = Nlx2MatVT(file,fieldSelection,extractHeader,extractMode);
                ind = find(x == 0);
                
                t(ind) = [];
                x(ind) = [];
                y(ind) = [];
                
                % Adjust input data
                D = 95;
                xscale = D/(max(x)-min(x));
                yscale = D/(max(y)-min(y));
                x = x * xscale;
                y = y * yscale;
                
                % Convert timestamps to seconds
                t = t/1000000;
                % Do median filtering to supress off-values
                %[x, y, t] = medianFilter(x, y, t);
                % Smoothing the position samples with a simple mean filter
                for cc=8:length(x)-7
                    x(cc) = nanmean(x(cc-7:cc+7));
                    y(cc) = nanmean(y(cc-7:cc+7));
                end
                
                % With many off-values in a row the median filter doesn't remove all, must
                % apply an additional off-value filter.
                %[x, y, t] = offValueFilter(x, y, t);
                
                [x,y] = centreBox(x,y);
                x_cen = (max(x) - min(x))/2 + min(x);
                y_cen = (max(y) - min(y))/2 + min(y);
                x = x - x_cen;
                y = y - y_cen;
                
                %compute velocity
                vel = zeros(1,length(x));
                for v = 2:1:(length(x)-1)
                    vel(v) = sqrt((x(v+1)-x(v-1))^2 + (y(v+1)-y(v-1))^2)/(t(v+1)-t(v-1));
                end
                vel(end) = sqrt((x(end)-x(end-1))^2 + (y(end)-y(end-1))^2)/(t(length(x))-t(length(x)-1));
                vel(vel>=40) = 0.5*(vel(circshift((vel>=40),-3))+vel(circshift((vel>=40),3)));
                vel = smooth(vel,15);
            end
            
            efile = strcat(sessions{jj},'Events.nev');
            FS = [1 0 1 0 1]; EH = 0; EM = 1;
            [TimeStamps, TTLs, EventStrings] = Nlx2MatEV(efile,FS,EH,EM);
            TTLs(1) = []; TimeStamps(1) = []; TimeStamps = TimeStamps/1000000;
            if length(TTLs) == 0
                msgbox('No TTLs in input file!','ERROR');
            end
            
            
            filex = [sessions{jj},channels{ii}];
            [samples,ts,tt, Fs, bv, ir] = loadEEG2(filex);
            x = bv*samples;
            filey = [sessions{jj},channels{kk}];
            [samples,ts,tt, Fs, bv, ir] = loadEEG2(filey);
            y = bv*samples;
            %set recording channel against ground and then against average reference
            xcscnum = channels{jj}(4:end-4); % get  tt# corresponding to avgsjj
            xcscnum = str2num(xcscnum);
            xttnum = ceil(xcscnum/4);
            ycscnum = channels{kk}(4:end-4); % get  tt# corresponding to avgsjj
            ycscnum = str2num(ycscnum);
            yttnum = ceil(ycscnum/4);
            if xttnum == 6 || xttnum == 10
                xrfile = [sessions{jj},refs{xttnum}];
                [samples,ts,tt, Fs, bv, ir] = loadEEG2(xrfile);
                x = x - bv*samples;
            end
            if yttnum == 6 || yttnum == 10
                yrfile = [sessions{jj},refs{yttnum}];
                [samples,ts,tt, Fs, bv, ir] = loadEEG2(yrfile);
                y = y - bv*samples;
            end
            %         if ~strcmp(refs{xttnum},'G')
            %             xrfile = [sessions{jj},refs{xttnum}];
            %             [samples,ts,tt, Fs, bv, ir] = loadEEG2(xrfile);
            %             x = x + bv*samples;
            %         end
            %         if ~strcmp(refs{yttnum},'G')
            %             yrfile = [sessions{jj},refs{yttnum}];
            %             [samples,ts,tt, Fs, bv, ir] = loadEEG2(yrfile);
            %             y = y + bv*samples;
            %         end
            %     x = x - avref;
            %     y = y - avref;
            clear samples
            
            [ImX_tmp] = traces2ImX(x,y,freqVec,Fs,width);
            
            %Set all vars to common time
            ImX_tmp(:,tt<TimeStamps(1))=[];tt(tt<TimeStamps(1))=[];
            vel(t<TimeStamps(1)) = []; t(t<TimeStamps(1))=[];
            ImX_tmp(:,tt>TimeStamps(end))=[];tt(tt>TimeStamps(end))=[];
            vel(t>TimeStamps(end)) = [];t(t>TimeStamps(end))=[];
            tt = tt - min(TimeStamps);
            t = t - min(TimeStamps);
            TS = TimeStamps - min(TimeStamps);
            
            %resize ImX to one sample per vel
            ImX_add = zeros(size(ImX_tmp,1),length(vel));
            w = median(diff(t));
            for i = 1:length(t)
                ImX_add(:,i) = nanmean(ImX_tmp(:,abs(tt-t(i))<=w),2);
            end
            ImX_t = [ImX_t ImX_add];
            clear ImX_tmp ImX_add
            
            if kk == 2
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
            end
            
            for ts=2:length(TS)
                tidx = tt>TS(ts-1)&tt<=TS(ts);
                if TTLs(ts) == 0
                    if kk == 2
                        veloff = [veloff vel(tidx)];
                    end
                    ImXoff = [ImXoff ImX(:,tidx)];
                else
                    if kk == 2
                        velon = [velon vel(tidx)];
                    end
                    ImXon = [ImXon ImX(:,tidx)];
                end
            end
        end
        
        %compute velocity-matched wpli
        velbins = min([veloff velon]):max([veloff velon]);
        %     trat = round(length(veloff)/length(velon)); %control for time spent in on/off states
        for vb = 1:length(velbins)-1
            
            %cumpute number of velocity samples to take from each state
            voffidx = find(veloff >= velbins(vb) & veloff < velbins(vb+1));
            vonidx = find(velon >= velbins(vb) & velon < velbins(vb+1));
            minvels = min([length(voffidx) length(vonidx)]);
            %         minvels = min([length(voffidx) trat*length(vonidx)]);
            
            %pick random subsample from each state
            off_ss = randsample(length(voffidx),minvels);
            on_ss = randsample(length(vonidx),minvels);
            %         on_ss = randsample(length(vonidx),minvels/trat);
            if minvels >=1 %&& minvels/trat >=1
                ImXonss = [ImXonss ImXon(:,vonidx(on_ss))];
                ImXoffss = [ImXoffss ImXoff(:,voffidx(off_ss))];
            end
        end
        close all
        figure(1);
        outsum(:,1)   = nansum(ImXoffss,2);      % compute the sum;
        outsumW(:,1)  = nansum(abs(ImXoffss),2); % normalization of the WPLI
        %     outssq(:,1)   = nansum(ImXoffss.^2,2);
        outsum(:,2)   = nansum(ImXonss,2);      % compute the sum;
        outsumW(:,2)  = nansum(abs(ImXonss),2); % normalization of the WPLI
        %     outssq(:,2)   = nansum(ImXonss.^2,2);
        %     wpli = abs(outsum.^2 - outssq)./(outsumW.^2 - outssq); % do the pairwise thing in a handy way
        %     wpli = wpli.^(0.5);
        wpli = abs(outsum)./outsumW;
        plot(freqVec,wpli(:,1),'k');
        hold on
%         plot(freqVec,zeros(size(freqVec)),'--k');
        plot(freqVec,wpli(:,2),'g');
        xlabel('frequency (Hz)');ylabel('wpli');
        title(['Opto-triggered WPLI: ',channels{ii}(1:end-4),'-',channels{kk}(1:end-4)]);
        hold off
        
        %         data.freqVec = freqVec; data.light_off = wpli(:,1); data.light_on = wpli(:,2);   %save data as struct
        %         filename = sprintf('%s%s%s%s',savedir,strcat('\OptoWPLI','\'),channels{ii}(1:end-4),'_',channels{kk}(1:end-4),'.mat');
        %         save(filename,'-struct','data');
        figImage = sprintf('%s%s%s%s',savedir,strcat('\OptoWPLI','\'),channels{ii}(1:end-4),'_',channels{kk}(1:end-4),'.fig');
                bmpImage = sprintf('%s%s%s%s',savedir,strcat('\OptoWPLI','\'),channels{ii}(1:end-4),'_',channels{kk}(1:end-4),'.bmp');
        %         epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_'),channels{jj}(1:end-4),'.eps');
        saveas(gcf,figImage,'fig');
                saveas(gcf,bmpImage,'bmp');
        %         saveas(gcf,epsImage,'eps');
    end
end


%end

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




