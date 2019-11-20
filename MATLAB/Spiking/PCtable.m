function PCtable(mice)

%create directory to store files
    dd = cd; %data directory
    dirInfo = dir(cd);
    found = 0;
    for ss=1:size(dirInfo,1)
        if dirInfo(ss).isdir
            if strcmp(dirInfo(ss).name,strcat('PCtable','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(dd,strcat('\PCtable','\'));
    end

fid = fopen(mice,'r');
if fid == -1
    msgbox('Could not open the input file! Make sure the filname and path are correct.','ERROR');
end

% Get the mice to be analyzed from input file
fid = fopen(mice,'r');
i = 1;

while ~feof(fid)
    str = fgetl(fid);
    if ~strcmp(str(end),'\')
        str = strcat(str,'\');
    end
    Mice(i) = {str};
    i = i+1;
end
numMice = i-1;

%initialize table for single lap place field measurements
pctable = [];
pcnames = cell(1,23);
pcnames{1} = 'fsize';
pcnames{2} = 'numspks';
pcnames{3} = 'peakpos';
pcnames{4} = 'com';
pcnames{5} = 'firstspk';
pcnames{6} = 'lastspk';
pcnames{7} = 'lapwidth';
pcnames{8} = 'sd';
pcnames{9} = 'skewness';
pcnames{10} = 'FRAI';
pcnames{11} = 'median';
pcnames{12} = 'tenter';
pcnames{13} = 'texit';
pcnames{14} = 'velenter';
pcnames{15} = 'velexit';
pcnames{16} = 'numfields';
pcnames{17} = 'rewdist';
pcnames{18} = 'fskew';
pcnames{19} = 'time';
pcnames{20} = 'session';
pcnames{21} = 'day';
pcnames{22} = 'cell';
pcnames{23} = 'mouse';

%initialize table for spike waveform measurements
spktable = [];
spknames = cell(1,19);
spknames{1} = 'amp';
spknames{2} = 'zamp';
spknames{3} = 'namp';
spknames{4} = 'dpamp';
spknames{5} = 'preamp';
spknames{6} = 'prezamp';
spknames{7} = 'prenamp';
spknames{8} = 'predpamp';
spknames{9} = 'burstpos';
spknames{10} = 'burstlen';
spknames{11} = 'preisi';
spknames{12} = 'postisi';
spknames{13} = 'pos';
spknames{14} = 'vel';
spknames{15} = 'time';
spknames{16} = 'session';
spknames{17} = 'day';
spknames{18} = 'cell';
spknames{19} = 'mouse';
numhist = 50; %number of bins in spike history term
histbins=10.^linspace(log10(0.001+0.01),log10(4+0.01),numhist+1)-0.01;
Hist = [];

%initialize some other variables
cn = 0; %cell number
D = 96.5; %track diameter in cm
binsize = 1; %cm
rbinsize = binsize/(pi*D)*2*pi; %rad
numbins = round(2*pi/rbinsize);
grid = linspace(-pi+pi/numbins,pi-pi/numbins,numbins); %place field evaluation grid
kappa = 554; %place field smoothing factor (~5 cm FWHM)
minr = 1;minp = 0.1; %min rate and proportion of peak rate
pad = 0.2; %padding for place field boundaries (~5 cm)
vmin = 0; %min velocity (~5 cm/sec)
fieldSelection(1) = 1; % Timestamps
fieldSelection(2) = 0; % Extracted X
fieldSelection(3) = 0; % Extracted Y
fieldSelection(4) = 0; % Extracted Angel
fieldSelection(5) = 1; % Targets
fieldSelection(6) = 0; % Points
extractHeader = 0; % Do we return header 1 = Yes, 0 = No.
extractMode = 1; % Extract all data
FieldSelectionFlags=[1 1 1 1 1]; %for ntt file
HeaderExtractionFlag=1; %for ntt file
ExtractionModeVector=[]; %for ntt file

%loop through mice
for i = 1:numMice
    disp(sprintf('%s%s','Reading data for: ',Mice{i}));
    
    %create directory to store files
    dirInfo = dir(Mice{i});
    found = 0;
    for ss=1:size(dirInfo,1)
        if dirInfo(ss).isdir
            if strcmp(dirInfo(ss).name,strcat('PCtable','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(Mice{i},strcat('PCtable','\'));
    end
    
    % find and open the dates.txt file
    dates_txt = strcat(Mice{i},'dates.txt');
    fid = fopen(dates_txt,'r');
    j = 1;
    
    while ~feof(fid)
        str = fgetl(fid);
        if ~strcmp(str(end),'\')
            str = strcat(str,'\');
        end
        Dates(j) = {str};
        j = j+1;
    end
    numDates = j-1;
    
    %loop through days for mouse 'i'
    for j=1:numDates
        %set directory
        disp(Dates{j});
        newpath = strcat(Mice{i},Dates{j});
        cd(newpath);
        inFile = 'infile - Copy.txt';
%         inFile = 'inFile.txt';
        
        %read .txt files for csc's, tt's, and field nums
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
        
        % read the file names from the tt-file list
        ttid = fopen(ttList,'r');
        jj = 1;
        while ~feof(ttid)
            str = fgetl(ttid);
            cells(jj) = {str};
            jj = jj+1;
        end    
        if cells{1} == -1 %skip if no place cells on this day
           continue
        end
        numcells = jj-1;    
        
        Time = cell(numsessions,1);
        Pos = cell(numsessions,1);
        Vel = cell(numsessions,1);
        mt = zeros(numsessions,1); %start time for each session
        
        %get all the tracking info for each session
        for ii = 1:numsessions
            disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));
            % Get position data            
            file = strcat(sessions{ii},'vt1.nvt');
            [Time{ii}, targets] = Nlx2MatVT(file,fieldSelection,extractHeader,extractMode);
            % Convert timestamps to seconds
            Time{ii} = Time{ii}'/1000000;
            % Decode the target data
            [dTargets,tracking] = decodeTargets(targets);
            % Exctract position data from the target data
            [frontX,frontY,backX,backY] = extractPosition(dTargets,tracking);        
            % Set missing samples to NaN.
            frontX(frontX==0&frontY==0) = NaN;frontY(frontX==0&frontY==0) = NaN;            
            backX(backX==0&backY==0) = NaN;backY(backX==0&backY==0) = NaN;     
            %scale and center trajectory
            nb = 200;
            xymap = zeros(nb,nb);
            v = 0.05;
            mx = linspace(min(backX),max(backX),nb);
            my = linspace(min(backY),max(backY),nb);
            for xx = 1:length(mx) 
                dx = mx(xx)-backX;
                for yy = 1:length(my)                    
                    dy = my(yy)-backY;
                    xymap(xx,yy) = nansum(v^2/sqrt(2*pi)*exp(-0.5*(dx.*dx*v^2+dy.*dy*v^2)));
                end
            end
            regmax = imregionalmax(xymap,4);
            [xx,yy] = find(regmax);            
            [z, a, b, alpha] = fitellipse([mx(xx);my(yy)]); %'linear' to speed up
            %Translate
            frontX = frontX-z(1);backX = backX-z(1);
            frontY = frontY-z(2);backY = backY-z(2);
            %Rotate
            Q = [cos(-alpha), -sin(-alpha); sin(-alpha) cos(-alpha)];
            ftemp = Q*[frontX;frontY];
            btemp = Q*[backX;backY];
            %Scale
            ftemp(1,:) = 0.5*D/a*ftemp(1,:);
            ftemp(2,:) = 0.5*D/b*ftemp(2,:);
            btemp(1,:) = 0.5*D/a*btemp(1,:);
            btemp(2,:) = 0.5*D/b*btemp(2,:);
            %Rotate back to orginal orientation
            Q = [cos(alpha), -sin(alpha); sin(alpha) cos(alpha)];
            ftemp = Q*ftemp;
            btemp = Q*btemp;
            frontX = ftemp(1,:);frontY = ftemp(2,:);
            backX = btemp(1,:);backY = btemp(2,:);
            %Fill missing position samples using polar interpolation
            [fphase,fr] = interporPos(frontX,frontY,2,29.97);
            [bphase,br] = interporPos(backX,backY,2,29.97);
            indfb = isnan(fphase) & isnan(bphase);
            Time{ii}(indfb) = [];fphase(indfb) = [];bphase(indfb) = [];   
            Pos{ii} = circ_mean([fphase;bphase]);           
            %smooth and get derivatives            
            Pos{ii} = unwrap(Pos{ii});
            span = 6/(max(Time{ii})-min(Time{ii}));
            Pos{ii} = smooth(Pos{ii},span,'loess');
            Vel{ii} = abs(diff(Pos{ii}))./diff(Time{ii});
            Vel{ii} = [Vel{ii}(1);Vel{ii}];
            span = 6/(max(Time{ii})-min(Time{ii}));
            Vel{ii} = exp(smooth(log(Vel{ii}),span,'loess'));
            %remove nans and reshape
            Time{ii}(isnan(Vel{ii})) = [];Pos{ii}(isnan(Vel{ii})) = [];
            Vel{ii}(isnan(Vel{ii})) = [];
            Pos{ii} = reshape(wraptopi(Pos{ii}),length(Pos{ii}),1);
            Vel{ii} = reshape(Vel{ii},length(Vel{ii}),1);
            %make time start at 0 for each session
            mt(ii) = min(Time{ii}); Time{ii} = Time{ii}-mt(ii);
        end       
        
        %loop through sessions of day 'j'
        for ii = 1:numsessions
            %create image save folders
            temp = sessions{ii};
            dirInfo = dir(temp);
            found2 = 0;
            for kk=1:size(dirInfo,1)
                if dirInfo(kk).isdir
                    if strcmp(dirInfo(kk).name,strcat('PlaceFieldAnglesPR','\'))
                        found2 = 1;
                    end
                end
            end
            if found2==0
                mkdir(temp,strcat('PlaceFieldAnglesPR','\'));
            end            
            %reset cell number
            if ii ~= 1
               cn = cn - numcells; 
            end    
            
            %loop through cells
            for jj=1:numcells
                
                %update cell number
                cn = cn + 1;
                
                %load spike times/positions for current session
                disp(sprintf('%s%s','Cell: ',num2str(jj)));                
                tfile = [sessions{ii},cells{jj}];                
                [TS] = loadSpikes(tfile);
                TS = TS-mt(ii);
                [spk_phase] = spikePos(TS,Pos{ii},Time{ii}); 
                
                %get spike waveforms from ntt file
                tfileidx = find(tfile=='_',1,'last');
                nttfile = strcat(tfile(1:tfileidx-1),'.ntt');                
                [AllTS,~,~,AllWavePeaks,AllWaves,~]=Nlx2MatSpike(nttfile, FieldSelectionFlags,HeaderExtractionFlag,extractMode,ExtractionModeVector);
                AllTS=AllTS./1000000 - mt(ii);
                AllWavePeaks = AllWavePeaks(1:4,:); % the first four are the peaks on each channel

                %only keep spikes for this cell
                WavePeaks = zeros(4,length(TS));
%                 Waves = zeros(32,4,length(TS));
                for nts=1:length(TS)
                    [~,idx]=min(abs(TS(nts)-AllTS));
                    WavePeaks(:,nts)=AllWavePeaks(:,idx);
%                     Waves(:,:,nts)=AllWaves(:,:,idx);
                end
                clear AllWaves AllWavePeaks AllTS
%                 [~,idx] = max(mean(WavePeaks,2)); %channel with largest spikes
%                 WavePeaks = WavePeaks(idx,:);
%                 Waves = squeeze(Waves(:,idx,:));
                WavePeaks = sqrt(sum(WavePeaks.^2));

                %prepare candidate entries for spktable
                tmpspktable = zeros(length(WavePeaks),length(spknames));
                tmpHist = zeros(length(WavePeaks),numhist+1);
                tmpspktable(:,1) = WavePeaks;clear WavePeaks %amp
                tmpspktable(:,2) = zscore(tmpspktable(:,1)); %zamp
                isis = diff(TS); %interspike intervals
                firstidx = find(isis>0.01) + 1; %indices for the first spike in each burst
                firstidx = [1; firstidx];
                numspks = diff(firstidx); %number of spikes in each burst
                if firstidx(end) ~= length(TS)
                    numspks = [numspks; length(TS)-firstidx(end)+1]; 
                else
                    numspks = [numspks; 1];
                end
                tmpspktable(1,3) = 1; %namp assuming first spike is also first spike in burst 
                bp = 1; %indexes position within burst
                tmpspktable(1,9) = bp; %burstpos
                bidx = 1; %indexes current burst number
                tmpspktable(1,10) = numspks(bidx); %burstlen                
                for nts = 2:length(TS)
                   recidx = max(firstidx(firstidx<=nts)); %indexes most recent first spike in a burst
                   tmpspktable(nts,3) = tmpspktable(nts,1)/tmpspktable(recidx,1); %namp 
                   if sum(firstidx == nts) == 1 
                       bp = 1; %first spike in burst
                       bidx = bidx+1; %increment burst index
                       tmpspktable(nts,9) = bp; %burstpos
                       tmpspktable(nts,10) = numspks(bidx); %burstlen
                   else 
                       bp = bp+1; %not first spike in burst
                       tmpspktable(nts,9) = bp; %burstpos
                       tmpspktable(nts,10) = numspks(bidx); %burstlen
                   end
                   tshist = TS(nts) - TS;                   
                   tmpHist(nts,:) = histc(tshist,histbins);
                end
                
                tmpspktable(1,4) = nan;
                tmpspktable(2:end,4) = tmpspktable(2:end,1)./tmpspktable(1:end-1,1); %dpamp
                tmpspktable(1,5) = nan;
                tmpspktable(2:end,5) = tmpspktable(1:end-1,1); %preamp
                tmpspktable(1,6) = nan;
                tmpspktable(2:end,6) = tmpspktable(1:end-1,2); %prezamp
                tmpspktable(1,7) = nan;
                tmpspktable(2:end,7) = tmpspktable(1:end-1,3); %prenamp
                tmpspktable(1,8) = nan;
                tmpspktable(2:end,8) = tmpspktable(1:end-1,4); %predpamp
                tmpspktable(1,11) = nan;
                tmpspktable(2:end,11) = isis; %preisi
                tmpspktable(end,12) = nan;
                tmpspktable(1:end-1,12) = isis; %postisi
                for nts=1:length(TS)
                    [~,idx]=min(abs(TS(nts)-Time{ii}));
                    tmpspktable(nts,13) = Pos{ii}(idx); %pos
                    tmpspktable(nts,14) = Vel{ii}(idx); %vel
                end
                tmpspktable(:,15) = TS; %time
                tmpspktable(:,16) = ii; %session
                tmpspktable(:,17) = j; %day
                tmpspktable(:,18) = cn; %cell
                tmpspktable(:,19) = i; %mouse
                
                %load spike times/positions/velocities for all sessions combined              
                V = [];P = [];S = [];T = [];SV = [];
                for kk = 1:numsessions
                    V = [V; Vel{kk}];
                    P = [P Pos{kk}'];
                    T = [T; Time{kk}]; 
                    tfile = [sessions{kk},cells{jj}];                
                    [spkt] = loadSpikes(tfile);
                    spkt = spkt - mt(kk);
                    [spkp] = spikePos(spkt,Pos{kk},Time{kk});
                    S = [S;spkp];         
                    [spkv] = spikePos(spkt,Vel{kk},Time{kk});
                    SV = [SV;spkv];
                end
                
                %estimate reward location
                pbins = linspace(-pi,pi,360);
                vmap = zeros(length(pbins),1);
                for p=1:length(pbins)
                    vmap(p) = mean(V(abs(angle(exp(1j*P).*conj(exp(1j*pbins(p)*ones(size(P))))))<=1/360*pi));
                end
                [~,im] = min(vmap);
                rloc = pbins(im);
                
                %remove spikes within some distance of reward location to prevent detecting false field
%                 spk_temp = S(abs(angle(exp(1j*S).*conj(exp(1j*pbins(im)*ones(size(S))))))>0); %was 0.4
                
                %compute ratemap on remaining spikes
                [map,~,~] = circle_map(S,P,T',kappa,grid);              
%                 counts = hist(S,grid);
%                 occup = hist(P,grid);
%                 map = counts./(occup*median(diff(T)));
                
                %test plotting
%                 plot(grid,map);xlim([-pi,pi]);hold on
%                 plot(grid(plocs),map(plocs),'ob');
%                 plot([pbins(im) pbins(im)],[0 max(map)]);keyboard
%                 hold off
                
                %detect place fields
                [numfields, peaks, plocs] = detectfields(map);  
                if isempty(peaks)
                    continue
                end                
                
                %get field boundaries
                cidx = plocs(peaks==max(peaks));
                [minb, maxb] = fieldboundaries(map,grid,cidx,minp,minr,pad); 
                
                %min distance of field to reward location
                rinfield = 0;
                if minb < maxb
                    if rloc>=minb && rloc<=maxb
                        rinfield=1;
                    end
                else
                    if rloc>minb || rloc<maxb
                        rinfield=1;
                    end
                end
                if rinfield
                    rd = 0;
                else
                    rd = min([abs(angle(exp(1j*minb).*conj(exp(1j*rloc)))),...
                        abs(angle(exp(1j*maxb).*conj(exp(1j*rloc))))]);
                end
                
                %if too close and >1 field, use next biggest field
                if  rd < 0.1 && numfields > 1 && i > 2
                    plocs(peaks==max(peaks)) = [];
                    peaks(peaks==max(peaks)) = [];
                    cidx = plocs(peaks==max(peaks));
                    [minb, maxb] = fieldboundaries(map,grid,cidx,minp,minr,pad);
                end              
                
                %determine which spikes are inside the boundaries
                field_id = ones(length(spk_phase),1);
                fmap = map;rinfield = 0;
                if minb < maxb
                    field_id(spk_phase<minb) = 0;
                    field_id(spk_phase>maxb) = 0;
                else
                    field_id(spk_phase<minb & spk_phase>maxb) = 0;
                end
                
                %throw away spikes outside the boundaries
                field_spks = spk_phase(field_id == 1);
                field_ts = TS(field_id == 1);       
                tmpspktable = tmpspktable(field_id==1,:);
                tmpHist = tmpHist(field_id==1,:);
                tmpHist(:,2) = zscore(tmpHist(:,1));
                
                %plot and save field detection figs
                plot(Time{ii},Pos{ii},'-k'); hold on
                plot(TS,spk_phase,'.b','MarkerSize',10);
                plot(field_ts,field_spks,'.r','MarkerSize',10);
                plot([Time{ii}(1) Time{ii}(end)],[maxb maxb],'-b');
                plot([Time{ii}(1) Time{ii}(end)],[minb minb],'-b');
                if ~isempty(field_spks)
                    plot([Time{ii}(1) Time{ii}(end)],[circ_median(field_spks) circ_median(field_spks)],'-r');
                end
                plot([Time{ii}(1) Time{ii}(end)],[pbins(im) pbins(im)],'--b');hold off
                xlabel('time (sec)');ylabel('track angle');       
                xlim([0 600]);ylim([-pi pi])
                title(strcat('Cell ',num2str(jj),'-',cells{jj}))
                bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('PlaceFieldAnglesPR','\'),cells{jj}(1:end-2),'.bmp');
                f = getframe(gcf);
                [pic, ~] = frame2im(f);
                imwrite(pic,bmpImage,'bmp');

                %rotate everything to be centered on place field center
                field_mu = circ_mean(field_spks); %field center  
                R_mu = [cos(-field_mu),-sin(-field_mu);sin(-field_mu),cos(-field_mu)]; %rotation matrix
                field_spks = (R_mu*[cos(field_spks),sin(field_spks)]')';
                field_spks = atan2(field_spks(:,2),field_spks(:,1)); %rotated spikes used for analyses
                path = (R_mu*[cos(Pos{ii}');sin(Pos{ii}')])';
                path = atan2(path(:,2),path(:,1)); %rotated position           
                minb = (R_mu*[cos(minb);sin(minb)])';
                minb = atan2(minb(:,2),minb(:,1)); %rotated field boundary
                maxb = (R_mu*[cos(maxb);sin(maxb)])';
                maxb = atan2(maxb(:,2),maxb(:,1)); %rotated field boundary
                spk_phase = (R_mu*[cos(spk_phase),sin(spk_phase)]')';
                spk_phase = atan2(spk_phase(:,2),spk_phase(:,1)); %rotated spikes only used for test plotting
                inspks = field_spks; %rotated spikes only used for test plotting
                rpos = (R_mu*[cos(tmpspktable(:,13)');sin(tmpspktable(:,13)')])';
                tmpspktable(:,13) = atan2(rpos(:,2),rpos(:,1)); %rotates spike positions for spike amplitude table
                rloc = (R_mu*[cos(rloc),sin(rloc)]')';
                rloc = atan2(rloc(2),rloc(1)); %reward location
                
                %get field skew
                rgrid = (R_mu*[cos(grid);sin(grid)])';
                rgrid = atan2(rgrid(:,2),rgrid(:,1));
                [rgrid,ridx] = sort(rgrid);
                rmap = map(ridx);
                rmap = rmap(rgrid>=minb & rgrid <= maxb);
                rgrid = rgrid(rgrid>=minb & rgrid <= maxb)';
                ratesum = sum(rmap);
                mu = rgrid*rmap/ratesum;
                sd = sqrt((rgrid-mu).^2*rmap/sum(rmap));
                fskew = sum((rgrid-mu).^3.*rmap')/((ratesum)*sd^3);
%                 fskew = (1/ratesum)*sum((rmap'.*rgrid-mu).^3) / ...
%                     ((1/ratesum)*sum((rmap'.*rgrid-mu).^2))^(3/2); %field skewness
                
                %detect single laps
                [~,~,~,imin] = extrema(cos(path));
                imin = imin(cos(path(imin))<-0.95);
                pass_times = sort(Time{ii}(imin));
                pass_times(diff(pass_times)<15) = [];
                %plots for testing rotations and lap cutting algorithm
%                 plot(Time{ii},cos(path),'.k');hold on;%plot(Time{ii},bp,'b')
%                 plot(field_ts,cos(field_spks),'*g')
%                 plot(pass_times,-1*ones(size(pass_times)),'.r');hold off
%                 pause
                clear imin bp   
                
                %Determine if cell has enough passes that meet criteria
                % i.e., median velocity > ~5 cm/sec                
                num_passes = 0; %num passes that meet criteria
                pass_idx = zeros(size(pass_times)); %index of accepted passes         
                meanfirst = 0;meanlast = 0;np = 0;
                for p = 1:length(pass_times)+1
                    if p==1 %conditions for first lap
                        numsp = sum(field_ts<pass_times(p));
                        if numsp>1
                            mints = min(field_ts(field_ts<pass_times(p)));
                            maxts = max(field_ts(field_ts<pass_times(p)));
                            medianv = median(Vel{ii}(Time{ii}>=mints&Time{ii}<=maxts));
                            if medianv > vmin 
                                num_passes = num_passes + 1;
                                pass_idx(p) = p;
                            end
                            meanfirst = meanfirst + min(field_spks(field_ts<pass_times(p)));
                            meanlast = meanlast + max(field_spks(field_ts<pass_times(p)));np=np+1;
                        else
                            field_spks(field_ts<pass_times(p)) = nan;
                            field_ts(field_ts<pass_times(p)) = nan;
                        end
                    elseif p>1 && p<=length(pass_times) %conditions for middle laps
                        numsp = sum(field_ts>pass_times(p-1)&field_ts<pass_times(p));
                        if numsp>1
                            mints = min(field_ts(field_ts>pass_times(p-1)&field_ts<pass_times(p)));
                            maxts = max(field_ts(field_ts>pass_times(p-1)&field_ts<pass_times(p)));
                            medianv = median(Vel{ii}(Time{ii}>=mints&Time{ii}<=maxts));
                            if medianv > vmin 
                                num_passes = num_passes + 1;
                                pass_idx(p) = p;
                            end
                            meanfirst = meanfirst + min(field_spks(field_ts>pass_times(p-1)&field_ts<pass_times(p)));
                            meanlast = meanlast + max(field_spks(field_ts>pass_times(p-1)&field_ts<pass_times(p)));np=np+1;
                        else
                            field_spks(field_ts>pass_times(p-1)&field_ts<pass_times(p)) = nan;
                            field_ts(field_ts>pass_times(p-1)&field_ts<pass_times(p)) = nan;
                        end
                    elseif p>length(pass_times) %condtions for last lap    
                        numsp = sum(field_ts>pass_times(p-1));
                        if numsp>1
                            mints = min(field_ts(field_ts>pass_times(p-1)));
                            maxts = max(field_ts(field_ts>pass_times(p-1)));
                            medianv = median(Vel{ii}(Time{ii}>=mints&Time{ii}<=maxts));
                            if medianv > vmin 
                                num_passes = num_passes + 1;
                                pass_idx(p) = p;
                            end
                            meanfirst = meanfirst + min(field_spks(field_ts>pass_times(p-1)));
                            meanlast = meanlast + max(field_spks(field_ts>pass_times(p-1)));np=np+1;
                        else
                            field_spks(field_ts>pass_times(p-1)) = nan;
                            field_ts(field_ts>pass_times(p-1)) = nan;
                        end
                    end
                end
%                 meanfirst = meanfirst/np;
%                 meanlast = meanlast/np;
                
                %remove passes that didn't meet criteria
                field_spks(isnan(field_spks)) = [];
                field_ts(isnan(field_ts)) = [];
                pass_idx(pass_idx == 0) = [];
                
                %only include cells that have enough laps
                if num_passes > 3              
                    
                    %apply second rotation based on remaining spikes
                    field_mu = circ_median(field_spks);
                    field_var = max(field_spks)-min(field_spks);
                    R_mu = [cos(-field_mu),-sin(-field_mu);sin(-field_mu),cos(-field_mu)];
                    field_spks = (R_mu*[cos(field_spks),sin(field_spks)]')';
                    field_spks = atan2(field_spks(:,2),field_spks(:,1));
                    path = (R_mu*[cos(path),sin(path)]')';
                    path = atan2(path(:,2),path(:,1));                    
                    minb = (R_mu*[cos(minb);sin(minb)])';
                    minb = atan2(minb(:,2),minb(:,1));
                    maxb = (R_mu*[cos(maxb);sin(maxb)])';
                    maxb = atan2(maxb(:,2),maxb(:,1));
                    spk_phase = (R_mu*[cos(spk_phase),sin(spk_phase)]')';
                    spk_phase = atan2(spk_phase(:,2),spk_phase(:,1));
                    inspks = (R_mu*[cos(inspks),sin(inspks)]')';
                    inspks = atan2(inspks(:,2),inspks(:,1));
                    rpos = (R_mu*[cos(tmpspktable(:,13)');sin(tmpspktable(:,13)')])';
                    tmpspktable(:,13) = atan2(rpos(:,2),rpos(:,1)); %rotates spike positions for spike amplitude table
                    rloc = (R_mu*[cos(rloc),sin(rloc)]')';
                    rloc = atan2(rloc(2),rloc(1)); %reward location
                    
                    %add entries to spktable                    
                    tmpspktable(:,13) = zscore(tmpspktable(:,13));
                    spktable = [spktable; tmpspktable]; 
                    Hist = [Hist; tmpHist];
                    clear tmpspktable tmpHist
                    
%                     %plot single cell for SfN poster (or general test plotting)
%                     sp = 60; %load MyColormaps; 
%                     cm = colormap; 
%              
%                     if ii == 1
%                         tm = max(Time{ii});
%                     elseif ii == 2
%                         tm = max(Time{ii-1})+max(Time{ii})+sp;
%                     else
%                         tm = max(Time{ii-2})+max(Time{ii-1})+max(Time{ii})+2*sp;
%                     end
%                     plot(Time{ii}+tm-max(Time{ii}),path,'.k','MarkerSize',0.5);hold on
%                     plot([tm+0.5*sp tm+0.5*sp],[-pi pi],'-k');
%                     plot(Time{ii}+tm-max(Time{ii}),zeros(size(Time{ii})),'--k');
%                     plot(Time{ii}+tm-max(Time{ii}),maxb*ones(size(Time{ii})),'-k');
%                     plot(Time{ii}+tm-max(Time{ii}),minb*ones(size(Time{ii})),'-k');
%                     plot(ints+tm-max(Time{ii}),inspks,'ok','MarkerSize',3,'MarkerFaceColor',[0.5 0.5 0.5]);
%                     ylim([-pi pi])
%                     xlim([0 max(Time{ii}+tm-max(Time{ii}))])

                    %loop through accepted passes
                    for p = pass_idx'
                        %conditions for first lap
                        if p==1 && sum(field_ts<pass_times(p)) > 0 && nanmean(cos(path(1:20))) < 0.5 ...
                                && sum(diff(Time{ii}(path>min(field_spks)&path<0&Time{ii}<pass_times(p)))) > 0 &&...
                                sum(diff(Time{ii}(path<max(field_spks)&path>0&Time{ii}<pass_times(p)))) > 0                      
                            
                            %get single pass ratemap
                            [passmap,~,~] = circle_map(field_spks(field_ts<pass_times(p)),path(Time{ii}<pass_times(p))',Time{ii}(Time{ii}<pass_times(p))',kappa,grid);
%                             counts = hist(field_spks(field_ts<pass_times(p)),grid);
%                             occup = hist(path(Time{ii}<pass_times(p))',grid);
%                             passmap = (counts./(occup*median(diff(T))))';
%                             plot(grid,passmap);xlim([-pi/2 pi/2]);pause()
                            
                            %add entries to pctable
                            pctable = [pctable; zeros(1,length(pcnames))];
                            pctable(end,1) = sum(passmap); %field size
                            pctable(end,2) = sum(field_ts<pass_times(p)); %numspks
                            if pctable(end,2) < 1 %need at least one spike for these measures
                                pctable(end,3) = nan; %peakpos                            
                                pctable(end,4) = nan; %com
                                pctable(end,11) = nan; %median
                            else
                                pctable(end,3) = grid(passmap==max(passmap)); %peakpos                            
                                pctable(end,4) = circ_mean(grid',passmap); %com
                                medidx = find(cumsum(passmap)/sum(passmap)<=0.5,1,'last');
                                pctable(end,11) = grid(medidx); %median
                            end                            
                            if pctable(end,2) < 2 %need at least two spikes for these measures
                                pctable(end,5) = nan; %firstspk
                                pctable(end,6) = nan; %lastspk
                                pctable(end,7) = nan; %lapwidth                                
                            else
                                pctable(end,5) = min(field_spks(field_ts<pass_times(p)));%-meanfirst; %firstspk
                                pctable(end,6) = max(field_spks(field_ts<pass_times(p)));%-meanlast; %lastspk
                                pctable(end,7) =  max(field_spks(field_ts<pass_times(p)))- min(field_spks(field_ts<pass_times(p))); %lapwidth                                
                            end 
                            if pctable(end,2) < 3 %need at least three spikes for these measures
                                pctable(end,8) = nan; %sd
                                pctable(end,9) = nan; %skewness
                                pctable(end,10) = nan; %FRAI
                            else
                                %compute skewness
                                ratesum = sum(rmap);
                                mu = rgrid*rmap/ratesum;
                                sd = sqrt((rgrid-mu).^2*rmap/sum(rmap));
                                fskew = sum((rgrid-mu).^3.*rmap')/((ratesum)*sd^3);
                                
                                ratesum = sum(passmap);
                                mu = grid*passmap/ratesum;
                                sd = sqrt((grid-mu).^2*passmap/ratesum);
                                pctable(end,8) = sd; %sd
                                pctable(end,9) = sum((grid-mu).^3.*passmap')/((ratesum)*sd^3); %skewness
%                                 pctable(end,8) = ((1/ratesum)*(sum((passmap'.*grid-mu).^3))) / ...
%                                     ((1/ratesum)*sum((passmap'.*grid-mu).^2))^(3/2); %skewness
                                %compute FRAI
%                                 isis = diff(field_ts(field_ts<pass_times(p))); %interspike intervals for pass
%                                 fone = 1/mean(isis(1:floor(length(isis)/2))); %firing rate in first half of spiking
%                                 ftwo = 1/mean(isis(ceil(length(isis)/2)+1:end)); %firing rate in second half of spiking
                                lapspks = field_ts(field_ts<pass_times(p));
                                medt = median(lapspks);
                                fone = sum(lapspks<=medt)/(medt-min(lapspks));
                                ftwo = sum(lapspks>=medt)/(max(lapspks)-medt);
                                pctable(end,10) = (fone-ftwo)/(fone+ftwo); %FRAI
                            end
                            pctable(end,12) = sum(diff(Time{ii}(path>min(field_spks)&path<0&Time{ii}<pass_times(p)))); %tenter
                            pctable(end,13) = sum(diff(Time{ii}(path<max(field_spks)&path>0&Time{ii}<pass_times(p)))); %texit
                            pctable(end,14) = abs(min(field_spks))/sum(diff(Time{ii}(path>min(field_spks)&path<0&Time{ii}<pass_times(p)))); %velenter
                            pctable(end,15) = max(field_spks)/sum(diff(Time{ii}(path<max(field_spks)&path>0&Time{ii}<pass_times(p)))); %velexit
                            pctable(end,16) = numfields; %numfields
                            pctable(end,17) = min(abs(angle(exp(1j*field_spks)*conj(exp(1j*rloc))))); %distance to reward
                            pctable(end,18) = fskew; %fskew
                            pctable(end,19) = median(field_ts(field_ts<pass_times(p))); %time
                            pctable(end,20) = ii; %session
                            pctable(end,21) = j; %day
                            pctable(end,22) = cn; %cell
                            pctable(end,23) = i; %mouse                                                       
                            
                            %test plotting
%                             cmidx = round(crrats(end,1)*256);
%                             if cmidx == 0
%                                 cmidx = 1;
%                             elseif cmidx > 256
%                                 cmidx = 256;
%                             end
%                             plot(field_ts(field_ts<pass_times(p))+tm-max(Time{ii}),field_spks(field_ts<pass_times(p)),'o','MarkerFaceColor',cm(cmidx,:),'MarkerEdgeColor','k')

                        %conditions for middle laps
                        elseif p>1 && p<=length(pass_times) && sum(field_ts>pass_times(p-1)&field_ts<pass_times(p)) > 0 ...
                                && sum(diff(Time{ii}(path>min(field_spks)&path<0&Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)))) > 0 ...
                            && sum(diff(Time{ii}(path<max(field_spks)&path>0&Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)))) > 0
                            
                            %get single pass ratemap
                            [passmap,~,~] = circle_map(field_spks(field_ts>pass_times(p-1)&field_ts<pass_times(p)),path(Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p))',...
                                Time{ii}(Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p))',kappa,grid);
%                             counts = hist(field_spks(field_ts>pass_times(p-1)&field_ts<pass_times(p)),grid);
%                             occup = hist(path(Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p))',grid);
%                             passmap = (counts./(occup*median(diff(T))))';keyboard
%                             passmap(isnan(passmap)) = 0;
%                             plot(grid,passmap);xlim([-pi/2 pi/2]);pause()                            
                            %add entries to pctable
                            pctable = [pctable; zeros(1,length(pcnames))];
                            pctable(end,1) = sum(passmap); %field size
                            pctable(end,2) = sum(field_ts>pass_times(p-1)&field_ts<pass_times(p)); %numspks
                            if pctable(end,2) < 1 %need at least one spike for these measures
                                pctable(end,3) = nan; %peakpos                            
                                pctable(end,4) = nan; %com
                                pctable(end,11) = nan; %median
                            else
                                pctable(end,3) = grid(passmap==max(passmap)); %peakpos                            
                                pctable(end,4) = circ_mean(grid',passmap); %com
                                medidx = find(cumsum(passmap)/sum(passmap)<=0.5,1,'last');
                                pctable(end,11) = grid(medidx); %median
                            end                            
                            if pctable(end,2) < 2 %need at least two spikes for these measures
                                pctable(end,5) = nan; %firstspk
                                pctable(end,6) = nan; %lastspk
                                pctable(end,7) = nan; %lapwidth                                
                            else
                                pctable(end,5) = min(field_spks(field_ts>pass_times(p-1)&field_ts<pass_times(p)));%-meanfirst; %firstspk
                                pctable(end,6) = max(field_spks(field_ts>pass_times(p-1)&field_ts<pass_times(p)));%-meanlast; %lastspk
                                pctable(end,7) =  max(field_spks(field_ts>pass_times(p-1)&field_ts<pass_times(p)))-min(field_spks(field_ts>pass_times(p-1)&field_ts<pass_times(p))); %lapwidth                                
                            end 
                            if pctable(end,2) < 3 %need at least three spikes for these measures
                                pctable(end,8) = nan; %sd
                                pctable(end,9) = nan; %skewness
                                pctable(end,10) = nan; %FRAI
                            else
                                %compute skewness
                                ratesum = sum(passmap);
                                mu = grid*passmap/ratesum;
                                sd = sqrt((grid-mu).^2*passmap/ratesum);
                                pctable(end,8) = sd; %sd 
                                pctable(end,9) = sum((grid-mu).^3.*passmap')/((ratesum)*sd^3); %skewness
%                                 pctable(end,8) = ((1/ratesum)*(sum((passmap'.*grid-mu).^3))) / ...
%                                     ((1/ratesum)*sum((passmap'.*grid-mu).^2))^(3/2); %skewness
                                %compute FRAI
%                                 isis = diff(field_ts(field_ts>pass_times(p-1)&field_ts<pass_times(p))); %interspike intervals for pass
%                                 fone = 1/mean(isis(1:floor(length(isis)/2))); %firing rate in first half of spiking
%                                 ftwo = 1/mean(isis(ceil(length(isis)/2)+1:end)); %firing rate in second half of spiking
                                lapspks = field_ts(field_ts>pass_times(p-1)&field_ts<pass_times(p));
                                medt = median(lapspks);
                                fone = sum(lapspks<=medt)/(medt-min(lapspks));
                                ftwo = sum(lapspks>=medt)/(max(lapspks)-medt);
                                pctable(end,10) = (fone-ftwo)/(fone+ftwo); %FRAI
                            end
                            pctable(end,12) = sum(diff(Time{ii}(path>min(field_spks)&path<0&Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)))); %tenter
                            pctable(end,13) = sum(diff(Time{ii}(path<max(field_spks)&path>0&Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)))); %texit
                            pctable(end,14) = abs(min(field_spks))/sum(diff(Time{ii}(path>min(field_spks)&path<0&Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)))); %velenter
                            pctable(end,15) = max(field_spks)/sum(diff(Time{ii}(path<max(field_spks)&path>0&Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)))); %velexit
                            pctable(end,16) = numfields; %numfields
                            pctable(end,17) = min(abs(angle(exp(1j*field_spks)*conj(exp(1j*rloc))))); %distance to reward
                            pctable(end,18) = fskew; %fskew
                            pctable(end,19) = median(field_ts(field_ts>pass_times(p-1)&field_ts<pass_times(p))); %time
                            pctable(end,20) = ii; %session
                            pctable(end,21) = j; %day
                            pctable(end,22) = cn; %cell
                            pctable(end,23) = i; %mouse                            
                            
                            %test plotting
%                             cmidx = round(crrats(end,1)*256);
%                             if cmidx == 0
%                                 cmidx = 1;
%                             elseif cmidx > 256
%                                 cmidx = 256;
%                             end
%                             plot(field_ts(field_ts>pass_times(p-1)&field_ts<pass_times(p))+tm-max(Time{ii}),field_spks(field_ts>pass_times(p-1)&field_ts<pass_times(p)),'o','MarkerFaceColor',cm(cmidx,:),'MarkerEdgeColor','k');                           

                        %conditions for last lap
                        elseif p>length(pass_times) && sum(field_ts>pass_times(p-1)) > 0 && nanmean(cos(path(end-20:end))) < 0.5 ...
                                && sum(diff(Time{ii}(path>min(field_spks)&path<0&Time{ii}>pass_times(p-1)))) > 0 && sum(diff(Time{ii}(path<max(field_spks)&path>0&Time{ii}>pass_times(p-1)))) > 0
                            
                           %get single pass ratemap
                            [passmap,~,~] = circle_map(field_spks(field_ts>pass_times(p-1)),path(Time{ii}>pass_times(p-1))',Time{ii}(Time{ii}>pass_times(p-1))',kappa,grid);
%                             counts = hist(field_spks(field_ts>pass_times(p-1)),grid);
%                             occup = hist(path(Time{ii}>pass_times(p-1))',grid);
%                             passmap = (counts./(occup*median(diff(T))))';
%                             plot(grid,passmap);xlim([-pi/2 pi/2]);pause()
                            
                            %add entries to pctable
                            pctable = [pctable; zeros(1,length(pcnames))];
                            pctable(end,1) = sum(passmap); %field size
                            pctable(end,2) = sum(field_ts>pass_times(p-1)); %numspks
                            if pctable(end,2) < 1 %need at least one spike for these measures
                                pctable(end,3) = nan; %peakpos                            
                                pctable(end,4) = nan; %com
                                pctable(end,11) = nan; %median
                            else
                                pctable(end,3) = grid(passmap==max(passmap)); %peakpos                            
                                pctable(end,4) = circ_mean(grid',passmap); %com
                                medidx = find(cumsum(passmap)/sum(passmap)<=0.5,1,'last');
                                pctable(end,11) = grid(medidx); %median
                            end                            
                            if pctable(end,2) < 2 %need at least two spikes for these measures
                                pctable(end,5) = nan; %firstspk
                                pctable(end,6) = nan; %lastspk
                                pctable(end,7) = nan; %lapwidth                                
                            else
                                pctable(end,5) = min(field_spks(field_ts>pass_times(p-1)));%-meanfirst; %firstspk
                                pctable(end,6) = max(field_spks(field_ts>pass_times(p-1)));%-meanlast; %lastspk
                                pctable(end,7) =  max(field_spks(field_ts>pass_times(p-1)))-min(field_spks(field_ts>pass_times(p-1))); %lapwidth                                
                            end 
                            if pctable(end,2) < 3 %need at least three spikes for these measures
                                pctable(end,8) = nan; %sd
                                pctable(end,9) = nan; %skewness
                                pctable(end,10) = nan; %FRAI
                            else
                                %compute skewness
                                ratesum = sum(passmap);
                                mu = grid*passmap/ratesum;
                                sd = sqrt((grid-mu).^2*passmap/ratesum);
                                pctable(end,8) = sd; %sd
                                pctable(end,9) = sum((grid-mu).^3.*passmap')/((ratesum)*sd^3); %skewness
%                                 pctable(end,8) = ((1/ratesum)*(sum((passmap'.*grid-mu).^3))) / ...
%                                     ((1/ratesum)*sum((passmap'.*grid-mu).^2))^(3/2); %skewness
                                %compute FRAI
%                                 isis = diff(field_ts(field_ts>pass_times(p-1))); %interspike intervals for pass
%                                 fone = 1/mean(isis(1:floor(length(isis)/2))); %firing rate in first half of spiking
%                                 ftwo = 1/mean(isis(ceil(length(isis)/2)+1:end)); %firing rate in second half of spiking
                                lapspks = field_ts(field_ts>pass_times(p-1));
                                medt = median(lapspks);
                                fone = sum(lapspks<=medt)/(medt-min(lapspks));
                                ftwo = sum(lapspks>=medt)/(max(lapspks)-medt);
                                pctable(end,10) = (fone-ftwo)/(fone+ftwo); %FRAI
                            end
                            pctable(end,12) = sum(diff(Time{ii}(path>min(field_spks)&path<0&Time{ii}>pass_times(p-1)))); %tenter
                            pctable(end,13) = sum(diff(Time{ii}(path<max(field_spks)&path>0&Time{ii}>pass_times(p-1)))); %texit
                            pctable(end,14) = abs(min(field_spks))/sum(diff(Time{ii}(path>min(field_spks)&path<0&Time{ii}>pass_times(p-1)))); %velenter
                            pctable(end,15) = max(field_spks)/sum(diff(Time{ii}(path<max(field_spks)&path>0&Time{ii}>pass_times(p-1)))); %velexit
                            pctable(end,16) = numfields; %numfields
                            pctable(end,17) = min(abs(angle(exp(1j*field_spks)*conj(exp(1j*rloc))))); %distance to reward
                            pctable(end,18) = fskew; %fskew
                            pctable(end,19) = median(field_ts(field_ts>pass_times(p-1))); %time
                            pctable(end,20) = ii; %session
                            pctable(end,21) = j; %day
                            pctable(end,22) = cn; %cell
                            pctable(end,23) = i; %mouse                            
                           
                            %test plotting
%                             cmidx = round(crrats(end,1)*256);
%                             if cmidx == 0
%                                 cmidx = 1;
%                             elseif cmidx > 256
%                                 cmidx = 256;
%                             end
%                             plot(field_ts(field_ts>pass_times(p-1))+tm-max(Time{ii}),field_spks(field_ts>pass_times(p-1)),'o','MarkerFaceColor',cm(cmidx,:),'MarkerEdgeColor','k');
                        end
                    end                    
                end                   
            end                
        end
    end
end

%save pctable as .csv
pctable = array2table(pctable,'VariableNames',pcnames);
pcfilename = sprintf('%s%s%s%s',dd,strcat('\PCtable','\'),'PCtable.csv');
writetable(pctable,pcfilename);

%save spktable as .csv
spktable = array2table(spktable,'VariableNames',spknames);
spkfilename = sprintf('%s%s%s%s',dd,strcat('\PCtable','\'),'spktable.csv');
writetable(spktable,spkfilename);

%save Hist as .csv
Hist = Hist(:,1:end-1);
Hist = Hist./diff(histbins);
histtable = array2table(Hist);
histfilename = sprintf('%s%s%s%s',dd,strcat('\PCtable','\'),'histtable.csv');
writetable(histtable,histfilename);

%_______________________________________________________________________
%other functions
%_______________________________________________________________________

% Finds the position of the spikes
function [pos] = spikePos(ts,phase,post)
N = length(ts);
pos = zeros(N,1);

for ii = 1:N
    tdiff = (post-ts(ii)).^2;
    [~,ind] = min(tdiff);
    pos(ii) = phase(ind(1));
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
function [theta,r] = interporPos(x,y,timeTreshold,sampRate)

theta = unwrap(atan2(y,x));
r = sqrt(x.^2+y.^2);

% Turn off warning
% warning('off','MATLAB:divideByZero');

% Number of samples that corresponds to the time threshold.
sampTreshold = floor(timeTreshold * sampRate);

% number of samples
numSamp = length(x);
% Find the indexes to the missing samples
% temp1 = 1./x;
% indt1 = isinf(temp1);
ind = isnan(x);
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
    theta(1:count) = theta(count);
    r(1:count) = r(count);
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
    theta(count:numSamp) = theta(count);
    r(count:numSamp) = r(count);
end

if change
    % Recalculate the missing samples
%     temp1 = 1./r;
%     indt1 = isinf(temp1);
    ind = isnan(r);
    % Missing samples are where both x and y are equal to zero
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
        theta(start) = theta(start-1);
        r(start) = r(start-1);
    else
        if count < sampTreshold
            % Last good position before lack of tracking
            theta1 = theta(start-1);
            r1 = r(start-1);
            % Next good position after lack of tracking
            theta2 = theta(stop);
            r2 = r(stop);
            % Calculate the interpolated positions
            theta(start:stop) = interp1([1,2],[theta1,theta2],1:1/count:2);
            r(start:stop) = interp1([1,2],[r1,r2],1:1/count:2);           
            % Increment the counter (avoid estimating allready estimated
            % samples)
            ii = ii+count;
        else
            % To many samples missing in a row and they are left as missing
            ii = ii+count;
        end
    end
end
theta = wraptopi(theta);

function [numfields, peaks, plocs] = detectfields(map)

[peaks,plocs] = findpeaks(map);
                [troughs,tlocs] = findpeaks(-map);
                troughs = abs(troughs);
                plocs = plocs(peaks>=2.5);
                peaks = peaks(peaks>=2.5);
                if isempty(peaks)
                    numfields = [];
                    peaks = [];
                    plocs = [];
                    return
                end
                if min(troughs(tlocs<plocs(1) | tlocs>plocs(end))) < 0.2*min([peaks(1) peaks(end)])
                else
                    if peaks(1)<peaks(end)
                        peaks(1) = [];
                        plocs(1) = [];
                    else
                        peaks(end) = [];
                        plocs(end) = [];
                    end
                end
                pidx = 1;
                for p = 1:length(peaks)-1
                    if min(troughs(tlocs>plocs(pidx) & tlocs<plocs(pidx+1))) < 0.2*min([peaks(pidx) peaks(pidx+1)])
                        pidx = pidx+1;
                    else
                        if peaks(pidx)<peaks(pidx+1)
                            peaks(pidx) = [];
                            plocs(pidx) = [];
                        else
                            peaks(pidx+1) = [];
                            plocs(pidx+1) = [];
                        end
                    end
                end
                numfields = length(peaks);

function [minb, maxb] = fieldboundaries(map,grid,cidx,minp,minr,pad)

minb = nan; maxb = nan; %place field boundaries for dominant field

count = 0;
for b = 1:length(grid)
    if cidx - b > 0
        bb = cidx - b;
    else
        bb = length(grid) - b + cidx;
    end
    if map(bb)<minp*map(cidx) %&& map(bb)<minr
        count = count + 1;
        if count == 10
            minb = grid(bb);
            break;
        end
    end
end
count = 0;
for b = 1:length(grid)
    if cidx + b <= length(grid)
        bb = cidx + b;
    else
        bb = cidx + b - length(grid);
    end
    if map(bb)<minp*map(cidx) %&& map(bb)<minr
        count = count + 1;
        if count == 10
            maxb = grid(bb);
            break;
        end
    end
end

% [~,troughs] = findpeaks(-map);
% premin = find(grid(troughs)<minb,1,'last');
% if ~isempty(premin)
%     minb = grid(troughs(premin));
% else
%     premin = max(troughs);
%     minb = grid(premin);
% end
% postmin = find(grid(troughs)>maxb,1,'first');
% if ~isempty(postmin)
%     maxb = grid(troughs(postmin));
% else
%     postmin = min(troughs);
%     maxb = grid(postmin);
% end

% if minb > -pi+pad
%     minb = minb - pad;
% else
%     minb = minb - pad + 2*pi;
% end
% if maxb < pi-pad
%     maxb = maxb + pad;
% else
%     maxb = maxb + pad - 2*pi;
% end

