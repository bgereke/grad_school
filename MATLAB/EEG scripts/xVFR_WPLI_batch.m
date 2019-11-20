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

function [xVFR] = xVFR_WPLI_batch(inFile,freqVec)

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

% Check if subdir for storing images are present. If not, it is
% created
dirInfo = dir(cd);
found = 0;
% foundm = 0;
% foundz = 0;
% founds = 0;
for kk=1:size(dirInfo,1)
    if dirInfo(kk).isdir
        if strcmp(dirInfo(kk).name,strcat('xVFR_WPLI','\'))
            found = 1;
        end
%         if strcmp(dirInfo(kk).name,strcat('mVFR_WPLI','\'))
%             foundm = 1;
%         end
%         if strcmp(dirInfo(kk).name,strcat('zVFR_WPLI','\'))
%             foundz = 1;
%         end
%         if strcmp(dirInfo(kk).name,strcat('sVFR_WPLI','\'))
%             founds = 1;
%         end
    end
end
if found==0
    mkdir(cd,strcat('xVFR_WPLI','\'));
end
% if foundm==0
%     mkdir(cd,strcat('mVFR_WPLI','\'));
% end
% if foundz==0
%     mkdir(cd,strcat('zVFR_WPLI','\'));
% end
% if founds==0
%     mkdir(cd,strcat('sVFR_WPLI','\'));
% end

vel = [];
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

for jj = 1:numchannels-1
    for kk = jj+1:numchannels
        ImX = [];
        for ii = 1:numsessions
            disp(sprintf('%s%s','Reading data for session: ',sessions{ii},'channels_',channels{jj},channels{kk}));
            if kk == 2             
                
                file = strcat(sessions{ii},'vt1.nvt');
                [t_t{ii}, x_t, y_t] = Nlx2MatVT(file,fieldSelection,extractHeader,extractMode);
                
                ind = find(x_t == 0);
                t_t{ii}(ind) = [];
                x_t(ind) = [];
                y_t(ind) = [];
                % Do median filtering to supress off-values
                %[x, y, t] = medianFilter(x, y, t);
                % Smoothing the position samples with a simple mean filter
                for cc=8:length(x_t)-7
                    x_t(cc) = nanmean(x_t(cc-7:cc+7));
                    y_t(cc) = nanmean(y_t(cc-7:cc+7));
                end
                y_t = -y_t + max(y_t)-min(y_t) + 2*max(y_t); % reflects posisitons so they are consistent with orientation in recording room
                [hx,xc] = hist(x_t,50); [hy,yc] = hist(y_t,50);
                th = 50;
                xmin = xc(find(hx>th,1,'first'));
                xmax = xc(find(hx>th,1,'last'));
                ymin = yc(find(hy>th,1,'first'));
                ymax = yc(find(hy>th,1,'last'));
                
                % Adjust input data
                scalex = 100/(xmax-xmin);
                scaley = 100/(ymax-ymin);
                x_t = x_t * scalex;
                y_t = y_t * scaley;
                xmin = xmin * scalex;
                ymin = ymin * scaley;
                xmax = xmax * scalex;
                ymax = ymax * scaley;
                xcen = xmin+(xmax-xmin)/2;
                ycen = ymin+(ymax-ymin)/2;
                x_t = x_t - xcen;
                y_t = y_t - ycen;
                
                % Convert timestamps to seconds
                t_t{ii} = t_t{ii}/1000000;
                
                %compute velocity and velocity map
                vel_t{ii} = zeros(1,length(x_t));
                for v = 2:1:(length(x_t)-1)
                    vel_t{ii}(v) = sqrt((x_t(v+1)-x_t(v-1))^2 + (y_t(v+1)-y_t(v-1))^2)/(t_t{ii}(v+1)-t_t{ii}(v-1));
                end
                vel_t{ii}(end) = sqrt((x_t(end)-x_t(end-1))^2 + (y_t(end)-y_t(end-1))^2)/(t_t{ii}(length(x_t))-t_t{ii}(length(x_t)-1));
                vel_t{ii}(1) = vel_t{ii}(2);
                vel_t{ii}(vel_t{ii}>=40) = 0.5*(vel_t{ii}(circshift((vel_t{ii}>=40),[-3,0]))+vel_t{ii}(circshift((vel_t{ii}>=40),[3,0])));
                vel_t{ii} = smooth(vel_t{ii},15);
                vel = [vel vel_t{ii}'];        
            end
            % Load data from the .ncs files, make plots, and store them
            xfile = [sessions{ii},channels{jj}];
            [samples,~,~,Fs,bv,~] = loadEEG2(xfile);
            chx{ii} = bv*samples;
            %set recording channel against ground and then against average reference
            xcscnum = channels{jj}(4:end-4); % get  tt# corresponding to avgsjj
            xcscnum = str2num(xcscnum);
            xttnum = ceil(xcscnum/4);
            if ~strcmp(refs{xttnum},'G')
                xrfile = [sessions{ii},refs{xttnum}];
                [samples,~,~,Fs,bv,~] = loadEEG2(xrfile);
                chx{ii} = chx{ii} + bv*samples;
            end
            yfile = [sessions{ii},channels{kk}];
            [samples,~,tt,Fs,bv,~] = loadEEG2(yfile);
            chy = bv*samples;
            ycscnum = channels{kk}(4:end-4); % get  tt# corresponding to avgsjj
            ycscnum = str2num(ycscnum);
            yttnum = ceil(ycscnum/4);
            if ~strcmp(refs{yttnum},'G')
                yrfile = [sessions{ii},refs{yttnum}];
                [samples,~,tt,Fs,bv,~] = loadEEG2(yrfile);
                chy = chy + bv*samples;
            end
            [ImX_tt] = traces2ImX(chx{ii},chy,freqVec,Fs,6);
            %resize ImX to one sample per vel
            ImX_t = zeros(size(ImX_tt,1),length(vel_t{ii}));
            w = median(diff(t_t{ii}));
            for i = 1:length(t_t{ii})
                ImX_t(:,i) = nanmean(ImX_tt(:,abs(tt-t_t{ii}(i))<=w),2);
            end
            ImX = [ImX ImX_t];
            clear x_t y_t ImX_t ImX_tt samples
        end      
        V = vel;
        V(isnan(ImX(1,:))) = [];
        ImX(:,isnan(ImX(1,:))) = [];
        [vh,vc] = hist(V,20);
        vmax = vc(vh<0.025*max(vh));
        vmax = min(vmax(vmax>10));
        vmax = min([vmax 30]);
        vbins = 0:0.5:vmax;
        lags = -30*5:30*5; %indices
        [xVFR,dxVFR,zxVFR,VFRsh] = xVFR_WPLI(ImX,V,vbins,lags);

        meanZ = zeros(length(freqVec),length(lags));
        meanSTD = zeros(length(freqVec),length(lags));
        for l = 1:length(lags)
           meanZ(:,l) = mean(abs(zxVFR(:,:,l)),2);
           meanSTD(:,l) = std(dxVFR(:,:,l),0,2);
        end
        figure(1)
        imagesc(lags/30,freqVec,meanZ,[0 max(max(meanZ(freqVec>50,:)))]);axis xy;colorbar;colormap hot;hold on;plot([0 0],[min(freqVec) max(freqVec)],'b')
        [mz, zidx] = max(meanZ,[],2);
        plot(lags(zidx(mz>0))/30,freqVec(mz>0),'.g');hold off
        xlabel('lag (sec)');ylabel('frequency (Hz)');
        title(strcat('VFR mean |z-score| from baseline: ',channels{jj}(1:end-4),'x',channels{kk}(1:end-4)));
        figImage = sprintf('%s%s%s%s',strcat('xVFR_WPLI','\zscore_'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.fig');
        bmpImage = sprintf('%s%s%s%s',strcat('xVFR_WPLI','\zscore_'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.bmp');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');
%         figure(2)
%         diffZ = meanZ(:,lags>0)-fliplr(meanZ(:,lags<0));
%         imagesc(lags(lags>0)/30,freqVec,diffZ,[-max(max(diffZ(freqVec>50,:))) max(max(diffZ(freqVec>50,:)))]);axis xy;colorbar;colormap jet;
%         xlabel('lag (sec)');ylabel('frequency (Hz)');
%         title(strcat('<|z-score|>_r_e_t_r_o - <|z-score|>_p_r_o: ',channels{jj}(1:end-4),'x',channels{kk}(1:end-4)));
%         figImage = sprintf('%s%s%s%s',strcat('xVFR_WPLI','\difference_'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.fig');
%         bmpImage = sprintf('%s%s%s%s',strcat('xVFR_WPLI','\difference_'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.bmp');
%         saveas(gcf,figImage,'fig');
%         saveas(gcf,bmpImage,'bmp');
        figure(3)
        imagesc(lags/30,freqVec,meanSTD,[0 max(max(meanSTD(freqVec>50,:)))]);axis xy;colorbar;colormap hot;hold on;plot([0 0],[min(freqVec) max(freqVec)],'b')
        [ms, sidx] = max(meanSTD,[],2);
        plot(lags(sidx(ms>0))/30,freqVec(ms>0),'.g');hold off
        xlabel('lag (sec)');ylabel('frequency (Hz)');
        title(strcat('VFR std from baseline: ',channels{jj}(1:end-4),'x',channels{kk}(1:end-4)));
        figImage = sprintf('%s%s%s%s',strcat('xVFR_WPLI','\STD_'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.fig');
        bmpImage = sprintf('%s%s%s%s',strcat('xVFR_WPLI','\STD_'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.bmp');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');

        data.freq = freqVec; data.vbins = vbins; data.xvfr = xVFR; data.dxvfr = dxVFR;   %save data as struct
        data.zxvfr = zxVFR; data.vfrsh = VFRsh; data.meanZ = meanZ;
        filename = sprintf('%s%s%s%s',strcat('xVFR_WPLI','\'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.mat');
        save(filename,'-struct','data');
        
    end
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



