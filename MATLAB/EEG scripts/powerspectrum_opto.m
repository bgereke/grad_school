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

function powerspectrum_opto(inFile,freqVec,width,numperms)  % width - number of cycles in wavelet (> 5 advisable); order - whitening

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

% % read the file names of the references for the channels to be used for
% % averaging
% refid = fopen(refList,'r');
% jj = 1;
% while ~feof(refid)
%        str = fgetl(refid);
%        refs(jj) = {str};
%        jj = jj+1;
% end

% % read the file names of the channels to be used for averaging
% avgid = fopen(ch4avg,'r');
% jj = 1;
% while ~feof(avgid)
%        str = fgetl(avgid);
%        avgs(jj) = {str};
%        jj = jj+1;
% end
% numrefs = jj-1;

%make directory to save figures
savedir = cd;
dirInfo = dir(cd);
found = 0;
for kk=1:size(dirInfo,1)
    if dirInfo(kk).isdir
        if strcmp(dirInfo(kk).name,strcat('OptoPowerSpectra','\'))
            found = 1;
        end
    end
end
if found==0
    mkdir(savedir,strcat('OptoPowerSpectra','\'));
end

veloff = [];
velon = [];
ons = [];
offs = [];

for ii = 1%:numchannels
    
    %declare main vars    
    TFRoff = [];
    TFRon = [];
    ps  = zeros(length(freqVec),2);
    normVec = zeros(2,1);

    disp(sprintf('%s%i',' CSC ',ii, ' of ',numchannels));
    
    for jj=1:numsessions
        
        % Get position data
        % Set the field selection for reading the video files. 1 = Add parameter, 0 = skip
        % parameter
        if ii == 1
            
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
%             % Do median filtering to supress off-values
%             %[x, y, t{jj}] = medianFilter(x, y, t{jj});
%             % Smoothing the position samples with a simple mean filter
%             for cc=8:length(x)-7
%                 x(cc) = nanmean(x(cc-7:cc+7));
%                 y(cc) = nanmean(y(cc-7:cc+7));
%             end
%             
%             % With many off-values in a row the median filter doesn't remove all, must
%             % apply an additional off-value filter.
%             %[x, y, t{jj}] = offValueFilter(x, y, t{jj});
%             
%             [x,y] = centreBox(x,y);
%             x_cen = (max(x) - min(x))/2 + min(x);
%             y_cen = (max(y) - min(y))/2 + min(y);
%             x = x - x_cen;
%             y = y - y_cen;
%             
%             %compute velocity
%             vel{jj} = zeros(1,length(x));
%             for v = 2:1:(length(x)-1)
%                 vel{jj}(v) = sqrt((x(v+1)-x(v-1))^2 + (y(v+1)-y(v-1))^2)/(t{jj}(v+1)-t{jj}(v-1));
%             end
%             vel{jj}(end) = sqrt((x(end)-x(end-1))^2 + (y(end)-y(end-1))^2)/(t{jj}(length(x))-t{jj}(length(x)-1));
%             vel{jj}(vel{jj}>=40) = 0.5*(vel{jj}(circshift((vel{jj}>=40),[-3,0]))+vel{jj}(circshift((vel{jj}>=40),[3,0])));
%             vel{jj} = smooth(vel{jj},15);
            
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
%                 velon = [velon vel{jj}(tidxon{jj})'];
%                 veloff = [veloff vel{jj}(tidxoff{jj})'];
            end            
        end
        
        if jj == 1 || jj == 2
            file = [sessions{jj},channels{ii}];
            [samples,~,tt, Fs, bv,~] = loadEEG2(file);
            ch_X = bv*samples;
            ch_X = ch_X(1:16:end);
            tt = tt(1:16:end);
            
%             xcscnum = channels{ii}(4:end-4); % get  tt# corresponding to avgsjj
%             xcscnum = str2num(xcscnum);
%             xttnum = ceil(xcscnum/4);
%             if xttnum == 6 || xttnum == 10
%                 xrfile = [sessions{jj},refs{xttnum}];
%                 [samples,~,~,~, bv,~] = loadEEG2(xrfile);
%                 ch_X = ch_X - bv*samples;
%             end
            %set recording channel against ground and then against average reference
            %         cscnum = channels{ii}(4:end-4); % get  tt# corresponding to avgsjj
            %         cscnum = str2num(cscnum);
            %         ttnum = ceil(cscnum/4);
            %         if ~strcmp(refs{ttnum},'G')
            %             rfile = [sessions{jj},refs{ttnum}];
            %             [samples,ts,tt, Fs, bv, ir] = loadEEG2(rfile);
            %             ch_X = ch_X + bv*samples;
            %         end
            %         ch_X = ch_X - avref;
            clear samples
            %         rfile = [sessions{jj},'R2.ncs'];
            %         [samples,ts,tt, Fs, bv, ir] = loadEEG2(rfile);
            %         Ref = bv*samples;
            %         ch_X = ch_X + Ref;
            %         ch_Xw = prewhitening(ch_X,3);   % prewhiten signal for better visualisation
            TFR_tmp = traces2TFR_rev(ch_X,freqVec,Fs,width);
              
            %resize TFR to one sample per vel
            TFR_t = zeros(size(TFR_tmp,1),length(t{jj}));
            tt = tt - mt{jj};
            w = median(diff(t{jj}));
            for i = 1:length(t{jj})
                TFR_t(:,i) = nanmean(TFR_tmp(:,abs(tt-t{jj}(i))<=w),2);
            end
            clear TFR_tmp

            TFRon = [TFRon TFR_t(:,tidxon{jj})];
            TFRoff = [TFRoff TFR_t(:,tidxoff{jj})];
            
            
            efile = strcat(sessions{jj},'Events.nev');
            FS = [1 0 1 0 1]; EH = 0; EM = 1;
            [TimeStamps, TTLs, EventStrings] = Nlx2MatEV(efile,FS,EH,EM);
            TTLs(1) = []; TimeStamps(1) = []; TimeStamps = TimeStamps/1000000;
            TS = TimeStamps - mt{jj};            
            for ts=2:length(TS)
                tidx = find(t{jj}>TS(ts-1)&t{jj}<=TS(ts));
                if TTLs(ts) == 0
                    offs = [offs mean(TFR_t(:,tidx),2)];
                else
                    ons = [ons mean(TFR_t(:,tidx),2)];
                end
            end           
            
            clear TFR_t
        end
    end    
    
%     %compute velocity-matched power spectra
%     velbins = min([veloff velon]):max([veloff velon]);
%     minvels = zeros(1,length(velbins)-1);    
% %     trat = round(length(veloff)/length(velon)); %control for time spent in on/off states
%     for vb = 1:length(velbins)-1
%     
%         %cumpute number of velocity samples to take from each state
%         voffidx = find(veloff >= velbins(vb) & veloff < velbins(vb+1));
%         vonidx = find(velon >= velbins(vb) & velon < velbins(vb+1));        
% %         minvels = min([length(voffidx) trat*length(vonidx)]);
%         minvels(vb) = min([length(voffidx) length(vonidx)]);
%         
%         %pick random subsample from each state
%         off_ss = randsample(length(voffidx),minvels(vb));
% %         on_ss = randsample(length(vonidx),minvels/trat);   
%         on_ss = randsample(length(vonidx),minvels(vb));  
%         normVec(1) = normVec(1)+minvels(vb);
%         ps(:,1) = ps(:,1) + nansum(TFRoff(:,voffidx(off_ss)),2);
%         normVec(2) = normVec(2)+minvels(vb);
%         ps(:,2) = ps(:,2) + nansum(TFRon(:,vonidx(on_ss)),2);
%     end    

    ps(:,1) = nanmean(TFRoff,2);
    ps(:,2) = nanmean(TFRon,2);
    
%     for jj=3:numsessions       
%         
%         file = [sessions{jj},channels{ii}];
%         [samples,~,tt, Fs, bv,~] = loadEEG2(file);
%         ch_X = bv*samples;
% %         xcscnum = channels{ii}(4:end-4); % get  tt# corresponding to avgsjj
% %         xcscnum = str2num(xcscnum);
% %         xttnum = ceil(xcscnum/4);
% %         if xttnum == 6 || xttnum == 10
% %             xrfile = [sessions{jj},refs{xttnum}];
% %             [samples,~,~,~, bv,~] = loadEEG2(xrfile);
% %             ch_X = ch_X - bv*samples;
% %         end
%         %set recording channel against ground and then against average reference
%         %         cscnum = channels{ii}(4:end-4); % get  tt# corresponding to avgsjj
%         %         cscnum = str2num(cscnum);
%         %         ttnum = ceil(cscnum/4);
%         %         if ~strcmp(refs{ttnum},'G')
%         %             rfile = [sessions{jj},refs{ttnum}];
%         %             [samples,ts,tt, Fs, bv, ir] = loadEEG2(rfile);
%         %             ch_X = ch_X + bv*samples;
%         %         end
%         clear samples
%         [TFR_tmp,~,freqVec] = traces2TFR_rev(ch_X,freqVec,Fs,width);
%         clear ch_X
%         %resize TFR to one sample per vel
%         TFR_t{jj} = zeros(size(TFR_tmp,1),length(vel{jj}));
%         tt = tt - mt{jj};
%         w = median(diff(t{jj}));
%         for i = 1:length(t{jj})
%             TFR_t{jj}(:,i) = nanmean(TFR_tmp(:,abs(tt-t{jj}(i))<=w),2);
%         end
%         clear TFR_tmp
% 
%     end
    
%     firstlight = linspace(0,30,numperms);
%     p_onps = zeros(length(freqVec),numperms);
%     p_offps = zeros(length(freqVec),numperms);    
%     
%     for b = 1:numperms
%         TFRon = []; TFRoff = [];
%         velon = []; veloff = [];          
%         
%         for jj = 3:numsessions
%             tidxon{jj} = []; tidxoff{jj} = [];
%             TS = firstlight(b):30:max(t{jj});
%             pTTLs = ones(size(TS));
%             pTTLs(1:2:end) = 0;
%             
%             for ts=2:length(TS)
%                 tidx = find(t{jj}>TS(ts-1)&t{jj}<=TS(ts));
%                 if pTTLs(ts) == 0
%                     tidxoff{jj} = [tidxoff{jj} tidx];
%                 else
%                     tidxon{jj} = [tidxon{jj} tidx];
%                 end
%             end
%             
%             velon = [velon vel{jj}(tidxon{jj})'];
%             veloff = [veloff vel{jj}(tidxoff{jj})'];
%             TFRon = [TFRon TFR_t{jj}(:,tidxon{jj})];
%             TFRoff = [TFRoff TFR_t{jj}(:,tidxoff{jj})];
%             
%         end       
%         
% %         for vb = 1:length(velbins)-1  
% % 
% %             %cumpute number of velocity samples to take from each state
% %             voffidx = find(veloff >= velbins(vb) & veloff < velbins(vb+1));
% %             vonidx = find(velon >= velbins(vb) & velon < velbins(vb+1));
% %             
% % %             off_ss = []; on_ss = [];
% %             minvels(vb) = min([length(voffidx) length(vonidx)]);
% %             
% %             %pick random subsample from each state
% %             off_ss = randsample(voffidx,minvels(vb));
% %             on_ss = randsample(vonidx,minvels(vb));
% %             
% % %             if ~isempty(voffidx) && ~isempty(vonidx)
% % %                 if minvels(vb) > length(voffidx)
% % %                     off_ss = [off_ss randsample(voffidx,minvels(vb),true)];
% % %                 else
% % %                     off_ss = [off_ss randsample(voffidx,minvels(vb))];
% % %                 end
% % %                 if minvels(vb) > length(vonidx)
% % %                     on_ss = [on_ss randsample(vonidx,minvels(vb),true)];
% % %                 else
% % %                     on_ss = [on_ss randsample(vonidx,minvels(vb))];
% % %                 end
% % %             end
% %             
% %             %pick random subsample from each state
% %             p_offps(:,b) = p_offps(:,b) + nansum(TFRoff(:,off_ss),2);
% %             p_onps(:,b) = p_onps(:,b) + nansum(TFRon(:,on_ss),2);
% %         end
%         p_offps(:,b) = nanmean(TFRoff,2);
%         p_onps(:,b) = nanmean(TFRon,2);
%     end
% 
%     pstd = std([p_onps-p_offps p_offps-p_onps],0,2);
    
    close all
    figure(1); 
%     subplot(1,2,1)
%     [ci] = prctile([(p_onps-p_offps)./p_offps*100 (p_offps-p_onps)./p_onps*100]',[2.5 97.5]);
%     shadedplot(freqVec, ci(2,:), ci(1,:),[1 0 0],[1 0.25 0.25]); alpha(0.5);hold on
%     plot(freqVec,mean([p_onps-p_offps p_offps-p_onps],2),'-r')
    plot(freqVec,log(ps(:,1)),'-k');hold on
    plot(freqVec,log(ps(:,2)),'-g')
    xlabel('frequency (Hz)')
    ylabel('log(power)')
    title(strcat('CSC ',channels{ii}(6:end-4)));
    
%     subplot(1,2,2)
%     [ci] = prctile([p_onps-p_offps p_offps-p_onps]'./repmat(pstd,1,2*numperms)',[2.5 97.5]);
%     shadedplot(freqVec, ci(2,:), ci(1,:),[1 0 0],[1 0.25 0.25]); alpha(0.5);hold on
%     plot(freqVec,mean([p_onps-p_offps p_offps-p_onps]'./repmat(pstd,1,2*numperms)'),'-r')
%     plot(freqVec,(ps(:,2)-ps(:,1))./pstd,'-k')
%     xlabel('frequency (Hz)')
%     ylabel('z-score from control')
%     title(strcat('CSC ',channels{ii}(4:end-4)));
    
%     figImage = sprintf('%s%s%s%s',savedir,strcat('\OptoPowerSpectra','\'),channels{ii}(1:end-4),'.fig');
    bmpImage = sprintf('%s%s%s%s',savedir,strcat('\OptoPowerSpectra','\'),channels{ii}(1:end-4),'.bmp');
%     saveas(gcf,figImage,'fig');
    saveas(gcf,bmpImage,'bmp');
    
    data.ps = ps; data.freqVec = freqVec; data.ons = ons; data.offs = offs;   %save data as struct
    filename = sprintf('%s%s%s%s',savedir,strcat('\OptoPowerSpectra','\'),channels{ii}(1:end-4),'.mat');
    save(filename,'-struct','data');
    
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




