function Pro_Retro_normalized(mice,freqVec,width)

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

img_text = 'on';
D = 100;

%create directory to store files
    dd = cd; %data directory
    dirInfo = dir(cd);
    found = 0;
    for ss=1:size(dirInfo,1)
        if dirInfo(ss).isdir
            if strcmp(dirInfo(ss).name,strcat('Pro_Retro_normalized','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(dd,strcat('\Pro_Retro_normalized','\'));
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

%declare group vars
gpass_spks = []; %data,pass,cell,session,day,mouse
gpass_ts = []; %data,pass,cell,session,day,mouse
gpass_path = []; %data,pass,cell,session,day,mouse
gpass_t = []; %data,pass,cell,session,day,mouse
gpass_numspks = []; %data,pass,cell,session,day,mouse
gfield_vars = []; %data,cell,session,day,mouse
gcmodes = []; %data,pass,cell,session,day,mouse
gpass_ct = []; %data,pass,cell,session,day,mouse
gCdiffs = []; %data
gTdiffs = []; %data
gpass_vel = []; %data,pass,cell,session,day,mouse
gspk_vel = [];
gaPFR = [];
gpPFR = [];
grPFR = [];
gaweights = [];
gpweights = [];
grweights = [];

%declare mouse vars
numvbins = 50;
vmin = 3.5;
aVFR = zeros(length(freqVec),numvbins,numMice);
pVFR = zeros(length(freqVec),numvbins,numMice);
rVFR = zeros(length(freqVec),numvbins,numMice);
avweights = zeros(length(freqVec),numvbins,numMice); 
pvweights = zeros(length(freqVec),numvbins,numMice);
rvweights = zeros(length(freqVec),numvbins,numMice);
aspkWPLI = zeros(length(freqVec),numMice); 
pspkWPLI = zeros(length(freqVec),numMice);
rspkWPLI = zeros(length(freqVec),numMice);
aspkweights = zeros(length(freqVec),numMice);
pspkweights = zeros(length(freqVec),numMice);
rspkweights = zeros(length(freqVec),numMice);
midx = 1;

for i = 1:numMice
    disp(sprintf('%s%s','Reading data for: ',Mice{i}));
    
    %create directory to store files
    dirInfo = dir(Mice{i});
    found = 0;
    for ss=1:size(dirInfo,1)
        if dirInfo(ss).isdir
            if strcmp(dirInfo(ss).name,strcat('Pro_Retro_normalized','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(Mice{i},strcat('Pro_Retro_normalized','\'));
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
    
    %declars vars being processed by pass/session
    pass_spks = []; %data,pass,cell,session,day,mouse
    pass_ts = []; %data,pass,cell,session,day,mouse
    pass_path = []; %data,pass,cell,session,day,mouse
    pass_t = []; %data,pass,cell,session,day,mouse
    pass_numspks = []; %data,pass,cell,session,day,mouse
    field_vars = []; %data,cell,session,day,mouse
    cmodes = []; %data,pass,cell,session,day,mouse
    pass_ct = []; %data,pass,cell,session,day,mouse
    Cdiffs = []; %data
    Tdiffs = []; %data
    pass_vel = []; %data,pass,cell,session,day,mouse
    spk_vel = []; %data,pass,cell,session,day,mouse
    session_start_idx = 1;
    numbins = 100;    
    aPFR = zeros(length(freqVec),numbins,numDates*3); %assumes begin1-3
    pPFR = zeros(length(freqVec),numbins,numDates*3);
    rPFR = zeros(length(freqVec),numbins,numDates*3);
    aweights = zeros(length(freqVec),numbins,numDates*3); 
    pweights = zeros(length(freqVec),numbins,numDates*3);
    rweights = zeros(length(freqVec),numbins,numDates*3);
    sidx = 1;
    aspkImX = []; 
    rspkImX = [];
    pspkImX = [];
    aV = [];
    rV = [];
    pV = [];
    
    % Perform all analyses for each date
    for j=1:numDates
        %set directory
        disp(Dates{j});
        newpath = strcat(Mice{i},Dates{j});
        cd(newpath);
        inFile = 'infile.txt';
        
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
        numcells = jj-1;
        
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
        
        % read the file names from the csc-file list
        cscid = fopen(cscList,'r');
        kk = 1;
        
        while ~feof(cscid)
            str = fgetl(cscid);
            channels(kk) = {str};
            kk = kk+1;
        end
        numchannels = kk-1;
        
        %loop through sessions
        for ii = 1:numsessions
            disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));
            
            %declare vars being processed by session
            aImX = []; 
            aImX_phase = []; 
            rImX = []; 
            rImX_phase = []; 
            pImX = [];
            pImX_phase = [];             
            
            % Get position data
            
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
            %             y = -y + max(y)-min(y) + 2*max(y); % reflects posisitons so they are consistent with orientation in recording room
            [hx,xc] = hist(x,50); [hy,yc] = hist(y,50);
            th = 50;
            xmin = xc(find(hx>th,1,'first'));
            xmax = xc(find(hx>th,1,'last'));
            ymin = yc(find(hy>th,1,'first'));
            ymax = yc(find(hy>th,1,'last'));
            
            % Adjust input data
            scalex = D/(xmax-xmin);
            scaley = D/(ymax-ymin);
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
            phase = atan2(y,x);
            
            % Convert timestamps to seconds
            t = t'/1000000;
            
            %compute velocity
            xs = smooth(x,15);ys = smooth(y,15);
            vel = zeros(1,length(x));
            for v = 2:1:(length(x)-1)
                vel(v) = sqrt((xs(v+1)-xs(v-1))^2 + (ys(v+1)-ys(v-1))^2)/(t(v+1)-t(v-1));
            end
            vel(end) = sqrt((xs(end)-xs(end-1))^2 + (ys(end)-ys(end-1))^2)/(t(length(x))-t(length(x)-1));
            vel(vel>=40) = 0.5*(vel(circshift((vel>=40),[-3,0]))+vel(circshift((vel>=40),[3,0])));
            vel = smooth(vel,30);
            
            %load CSC's
            filex = [sessions{ii},channels{1}];
            [samples,ts,tt, Fs, bv, ir] = loadEEG2(filex);
            ch_X = bv*samples;
            filey = [sessions{ii},channels{2}];
            [samples,ts,tt, Fs, bv, ir] = loadEEG2(filey);
            ch_Y = bv*samples;           
            
            %set recording channel against ground
            if i < 2 % i.e. don't do this for mouse32
                xcscnum = channels{1}(4:end-4); % get  tt# corresponding to avgsjj
                xcscnum = str2num(xcscnum);
                xttnum = ceil(xcscnum/4);
                ycscnum = channels{2}(4:end-4); % get  tt# corresponding to avgsjj
                ycscnum = str2num(ycscnum);
                yttnum = ceil(ycscnum/4);
                if ~strcmp(refs{xttnum},'G')
                    xrfile = [sessions{ii},refs{xttnum}];
                    [samples,ts,tt, Fs, bv, ir] = loadEEG2(xrfile);
                    ch_X = ch_X + bv*samples;
                end
                if ~strcmp(refs{yttnum},'G')
                    yrfile = [sessions{ii},refs{yttnum}];
                    [samples,ts,tt, Fs, bv, ir] = loadEEG2(yrfile);
                    ch_Y = ch_Y + bv*samples;
                end
            end
            [ImX] = traces2ImX(ch_X,ch_Y,freqVec,Fs,width);
            signs = sign(mean(ImX,2));
            ImX = ImX.*repmat(signs,1,size(ImX,2));
            clear Ch_X Ch_Y samples x y signs
            [ImX_phase] = spikePos(tt,phase,t);
            
            %resize vel to one sample per tt            
            v = zeros(length(tt),1);
            [~, ttidx] = min((tt-t(1)).^2);
            v(1:ttidx) = vel(1);
            old = ttidx+1;
            for ti = 2:length(t)
                [~, ttidx] = min((tt-t(ti)).^2);
                v(old:ttidx) = linspace(vel(ti-1),vel(ti),ttidx-old+1);
                old = ttidx+1;
            end
            if old < length(tt)
                v(old:end) = vel(end);
            end
            
            %normalize times
            mt = min(tt); t = t-mt; tt = tt-mt; 
            
            %loop through cells
            for jj=1:numcells
                
                disp(sprintf('%s%s','Cell: ',num2str(jj)));                
                tfile = [sessions{ii},cells{jj}];                
                [TS] = loadSpikes(tfile);
                TS = TS-mt;
                [spk_phase] = spikePos(TS,phase,t); 
                [spkv] = spikePos(TS,vel,t); 
                spk_temp = spk_phase(spkv>5);
                
                %compute place field
                binsize = 5; %cm
                mapAxis = linspace(-pi,pi,pi*D/binsize);
                kappa = 100; % Smoothing factor when calculating the ratemap
                minr = 1; minp = 0.2; %min rate and proportion of peak rate
                [map,~,cidx] = circle_map(spk_temp,phase,t',kappa,mapAxis);
                %Set Boundaries ~Kevin's method
                minb = nan; maxb = nan;
                for b = 1:length(mapAxis)
                    if cidx - b > 0
                        bb = cidx - b;
                    else
                        bb = length(mapAxis) - b + cidx;
                    end
                    if map(bb)<minr && map(bb)<minp*map(cidx)
                        minb = mapAxis(bb);
                    end
                    if ~isnan(minb), break; end
                end
                for b = 1:length(mapAxis)
                    if cidx + b <= length(mapAxis)
                        bb = cidx + b;
                    else
                        bb = cidx + b - length(mapAxis);
                    end
                    if map(bb)<minr && map(bb)<minp*map(cidx)
                        maxb = mapAxis(bb);
                    end
                    if ~isnan(maxb), break; end
                end
                field_id = ones(length(spk_phase),1);
                if minb < maxb
                    field_id(spk_phase<minb) = 0;
                    field_id(spk_phase>maxb) = 0;
                else
                    field_id(spk_phase<minb & spk_phase>maxb) = 0;
                end            
                
                field_spks = spk_phase(field_id == 1);
                field_ts = TS(field_id == 1);
                
                %find field center and size
                field_mu = circ_median(field_spks);
                [field_var, ~] = circ_var(field_spks);
                
                %rotate spikes to be centered at 0
                R_mu = [cos(-field_mu),-sin(-field_mu);sin(-field_mu),cos(-field_mu)];
                field_spks = (R_mu*[cos(field_spks),sin(field_spks)]')';
                field_spks = atan2(field_spks(:,2),field_spks(:,1));
                
                %rotate animal trajectory accordingly and break into single passes
                path = (R_mu*[cos(phase);sin(phase)])';
                path = atan2(path(:,2),path(:,1));
                ImX_phase = (R_mu*[cos(ImX_phase),sin(ImX_phase)]')';
                ImX_phase = atan2(ImX_phase(:,2),ImX_phase(:,1));
%                 pass_times = t(abs(diff(path))>3*pi/2);
%                 bp = fftlowpass(cos(path),30,0.1,0.2);                 
                [~,~,~,imin] = extrema(cos(path));
                imin = imin(cos(path(imin))<-0.95);
                pass_times = sort(t(imin));               
                pass_times(diff(pass_times)<15) = [];
                %plots for testing rotations and lap cutting algorithm
%                 plot(t,cos(path),'k');hold on;%plot(t,bp,'b')
%                 plot(field_ts,cos(field_spks),'*g')
%                 plot(pass_times,-1*ones(size(pass_times)),'.r');hold off
%                 pause
                clear imin bp
                
                %Determine if cell has enough passes that meet criteria
                num_passes = 0; %num passes that meet criteria
                pass_idx = zeros(size(pass_times)); %index of accepted passes
                for p = 1:length(pass_times)+1
                    if p==1
                        numsp = sum(field_ts<pass_times(p+1));
                        if numsp>1
                            %                             minv = min(vel(t>min(field_ts(field_ts>pass_times(p-1)&field_ts<pass_times(p)))...
                            %                             &t<max(field_ts(field_ts>pass_times(p-1)&field_ts<pass_times(p)))));
                            minv = min(vel(t<pass_times(p+1)&path>min(field_spks)&path<max(field_spks)));
                            if minv>vmin
                                num_passes = num_passes + 1;
                                pass_idx(p) = p;
                            end
                        else
                            field_spks(field_ts<pass_times(p+1)) = nan;
                            field_ts(field_ts<pass_times(p+1)) = nan;
                        end
                    elseif p>1 && p<=length(pass_times)
                        numsp = sum(field_ts>pass_times(p-1)&field_ts<pass_times(p));
                        if numsp>1
                            %                             minv = min(vel(t>min(field_ts(field_ts>pass_times(p-1)&field_ts<pass_times(p)))...
                            %                             &t<max(field_ts(field_ts>pass_times(p-1)&field_ts<pass_times(p)))));
                            minv = min(vel(t>pass_times(p-1)&t<pass_times(p)&path>min(field_spks)&path<max(field_spks)));
                            if minv>vmin
                                num_passes = num_passes + 1;
                                pass_idx(p) = p;
                            end
                        else
                            field_spks(field_ts>pass_times(p-1)&field_ts<pass_times(p)) = nan;
                            field_ts(field_ts>pass_times(p-1)&field_ts<pass_times(p)) = nan;
                        end
                    elseif p==length(pass_times)    
                        numsp = sum(field_ts>pass_times(p));
                        if numsp>1
                            %                             minv = min(vel(t>min(field_ts(field_ts>pass_times(p-1)&field_ts<pass_times(p)))...
                            %                             &t<max(field_ts(field_ts>pass_times(p-1)&field_ts<pass_times(p)))));
                            minv = min(vel(t>pass_times(p)&path>min(field_spks)&path<max(field_spks)));
                            if minv>vmin
                                num_passes = num_passes + 1;
                                pass_idx(p) = p;
                            end
                        else
                            field_spks(field_ts>pass_times(p)) = nan;
                            field_ts(field_ts>pass_times(p)) = nan;
                        end
                    end
                end
                
                field_spks(isnan(field_spks)) = []; 
                field_ts(isnan(field_ts)) = [];
                pass_idx(pass_idx == 0) = [];
                
                if num_passes > 3
                    %apply second rotation based on remaining spikes
                    field_mu = circ_median(field_spks);
                    [field_var, ~] = circ_var(field_spks);
                    R_mu = [cos(-field_mu),-sin(-field_mu);sin(-field_mu),cos(-field_mu)];
                    field_spks = (R_mu*[cos(field_spks),sin(field_spks)]')';
                    field_spks = atan2(field_spks(:,2),field_spks(:,1));
                    ImX_phase = (R_mu*[cos(ImX_phase),sin(ImX_phase)]')';
                    ImX_phase = atan2(ImX_phase(:,2),ImX_phase(:,1));
                    path = (R_mu*[cos(path),sin(path)]')';
                    path = atan2(path(:,2),path(:,1));
                    %loop through accepted passes
                    for p = pass_idx'        
                        if p==1 && sum(field_ts<pass_times(p+1)) > 1 && abs(path(1)) > 2 && abs(path(find(t<pass_times(p+1),1,'last'))) > 2
                            
                            %get pass spikes
                            spks = field_spks(field_ts<pass_times(p+1));
                            
                            %determine coding mode
                            retro = 0; pro = 0;
                            if sum(spks>0)/length(spks) > 2/3
                                retro = 1;
                            elseif sum(spks<0)/length(spks) > 2/3
                                pro = 1;
                            end
                            
                            %get first spike, apply third rotation, and scale so that spikes fall between 0 and 1
                            fspk = min(spks);
                            R_mu = [cos(-fspk),-sin(-fspk);sin(-fspk),cos(-fspk)];
                            spks = (R_mu*[cos(spks),sin(spks)]')';
                            spks = atan2(spks(:,2),spks(:,1));
                            passImX_phase = (R_mu*[cos(ImX_phase),sin(ImX_phase)]')';
                            passImX_phase = atan2(passImX_phase(:,2),passImX_phase(:,1));
                            passpath = (R_mu*[cos(path),sin(path)]')';
                            passpath = atan2(passpath(:,2),passpath(:,1));
                            scale = 1/max(spks);
                            spks = scale*spks; passpath = scale*passpath;
                            passImX_phase = scale*passImX_phase;
                            
                            mint = t(find(t<pass_times(p+1),1,'first'));
                            maxt = t(find(t<pass_times(p+1),1,'last'));
                            
                            %process vars
                            pass_spks = [pass_spks; [spks p*ones(length(spks),1) jj*j*ones(length(spks),1) ii*ones(length(spks),1)...
                                j*ones(length(spks),1) i*ones(length(spks),1)]];
                            pass_ts = [pass_ts; [field_ts(field_ts>mint&field_ts<maxt) p*ones(sum(field_ts>mint&field_ts<maxt),1)...
                                jj*j*ones(sum(field_ts>mint&field_ts<maxt),1) ii*ones(sum(field_ts>mint&field_ts<maxt),1)...
                                j*ones(sum(field_ts>mint&field_ts<maxt),1) i*ones(sum(field_ts>mint&field_ts<maxt),1)]];
                            pass_ct = [pass_ct; [median(field_ts(field_ts>mint&field_ts<maxt)) p jj*j ii j i]];
                            pass_path = [pass_path; [passpath(t>mint&t<maxt) p*ones(sum(t>mint&t<maxt),1)...
                                jj*j*ones(sum(t>mint&t<maxt),1) ii*ones(sum(t>mint&t<maxt),1)...
                                j*ones(sum(t>mint&t<maxt),1) i*ones(sum(t>mint&t<maxt),1)]];
                            pass_t = [pass_t; [t(t>mint&t<maxt) p*ones(sum(t>mint&t<maxt),1)...
                                jj*j*ones(sum(t>mint&t<maxt),1) ii*ones(sum(t>mint&t<maxt),1)...
                                j*ones(sum(t>mint&t<maxt),1) i*ones(sum(t>mint&t<maxt),1)]];
                            pass_numspks = [pass_numspks; [length(spks) p jj*j ii j i]];    
                            field_vars = [field_vars; [field_var jj*j ii j i]];
                            pass_vel = [pass_vel;[vel(t>mint&t<maxt,1) p*ones(sum(t>mint&t<maxt),1) jj*j*ones(sum(t>mint&t<maxt),1)...
                                ii*ones(sum(t>mint&t<maxt),1) j*ones(sum(t>mint&t<maxt),1) i*ones(sum(t>mint&t<maxt),1)]];  
                            spk_vel = [spk_vel; [vel(t>min(field_ts(field_ts>mint))&t<max(field_ts(field_ts<maxt)),1)...
                                p*ones(sum(t>min(field_ts(field_ts>mint))&t<max(field_ts(field_ts<maxt))),1)...
                                jj*j*ones(sum(t>min(field_ts(field_ts>mint))&t<max(field_ts(field_ts<maxt))),1)...
                                ii*ones(sum(t>min(field_ts(field_ts>mint))&t<max(field_ts(field_ts<maxt))),1)...
                                j*ones(sum(t>min(field_ts(field_ts>mint))&t<max(field_ts(field_ts<maxt))),1)...
                                i*ones(sum(t>min(field_ts(field_ts>mint))&t<max(field_ts(field_ts<maxt))),1)]];                           

                            if retro == 1
                                cmodes = [cmodes;[1 p jj*j ii j i]];
                                rImX = [rImX ImX(:,tt>mint&tt<maxt)];
                                rImX_phase = [rImX_phase;passImX_phase(tt>mint&tt<maxt)];
                                rspkImX = [rspkImX ImX(:,tt>min(field_ts(field_ts>mint))&tt<max(field_ts(field_ts<maxt)))];
                                rV = [rV;v(tt>min(field_ts(field_ts>mint))&tt<max(field_ts(field_ts<maxt)))]; 
                            elseif pro == 1
                                cmodes = [cmodes;[-1 p jj*j ii j i]];
                                pImX = [pImX ImX(:,tt>mint&tt<maxt)];
                                pImX_phase = [pImX_phase;passImX_phase(tt>mint&tt<maxt)];
                                pspkImX = [pspkImX ImX(:,tt>min(field_ts(field_ts>mint))&tt<max(field_ts(field_ts<maxt)))];
                                pV = [pV;v(tt>min(field_ts(field_ts>mint))&tt<max(field_ts(field_ts<maxt)))]; 
                            else
                                cmodes = [cmodes;[nan p jj*j ii j i]];
                                aImX = [aImX ImX(:,tt>mint&tt<maxt)];
                                aImX_phase = [aImX_phase;passImX_phase(tt>mint&tt<maxt)];
                                aspkImX = [aspkImX ImX(:,tt>min(field_ts(field_ts>mint))&tt<max(field_ts(field_ts<maxt)))];                                    
                                aV = [aV;v(tt>min(field_ts(field_ts>mint))&tt<max(field_ts(field_ts<maxt)))]; 
                            end                                        
                            
                        elseif p>1 && p<=length(pass_times) && sum(field_ts>pass_times(p-1)&field_ts<pass_times(p)) > 1 %if it's a full lap and there are at least two spikes                             
                            
                            %get pass spikes
                            spks = field_spks(field_ts>pass_times(p-1)&field_ts<pass_times(p));
                            
                            %determine coding mode
                            retro = 0; pro = 0;
                            if sum(spks>0)/length(spks) > 2/3
                                retro = 1;
                            elseif sum(spks<0)/length(spks) > 2/3
                                pro = 1;
                            end
                            
                            %get first spike, apply third rotation, and scale so that spikes fall between 0 and 1
                            fspk = min(spks);
                            R_mu = [cos(-fspk),-sin(-fspk);sin(-fspk),cos(-fspk)];
                            spks = (R_mu*[cos(spks),sin(spks)]')';
                            spks = atan2(spks(:,2),spks(:,1));
                            passImX_phase = (R_mu*[cos(ImX_phase),sin(ImX_phase)]')';
                            passImX_phase = atan2(passImX_phase(:,2),passImX_phase(:,1));
                            passpath = (R_mu*[cos(path),sin(path)]')';
                            passpath = atan2(passpath(:,2),passpath(:,1));
                            scale = 1/max(spks);
                            spks = scale*spks; passpath = scale*passpath;
                            passImX_phase = scale*passImX_phase;
                            
                            mint = min(t(t>=pass_times(p-1)&t<=pass_times(p)));
                            maxt = max(t(t>=pass_times(p-1)&t<=pass_times(p)));
                            
                            %process vars
                            pass_spks = [pass_spks; [spks p*ones(length(spks),1) jj*j*ones(length(spks),1) ii*ones(length(spks),1)...
                                j*ones(length(spks),1) i*ones(length(spks),1)]];
                            pass_ts = [pass_ts; [field_ts(field_ts>mint&field_ts<maxt) p*ones(sum(field_ts>mint&field_ts<maxt),1)...
                                jj*j*ones(sum(field_ts>mint&field_ts<maxt),1) ii*ones(sum(field_ts>mint&field_ts<maxt),1)...
                                j*ones(sum(field_ts>mint&field_ts<maxt),1) i*ones(sum(field_ts>mint&field_ts<maxt),1)]];
                            pass_ct = [pass_ct; [median(field_ts(field_ts>mint&field_ts<maxt)) p jj*j ii j i]];
                            pass_path = [pass_path; [passpath(t>mint&t<maxt) p*ones(sum(t>mint&t<maxt),1)...
                                jj*j*ones(sum(t>mint&t<maxt),1) ii*ones(sum(t>mint&t<maxt),1)...
                                j*ones(sum(t>mint&t<maxt),1) i*ones(sum(t>mint&t<maxt),1)]];
                            pass_t = [pass_t; [t(t>mint&t<maxt) p*ones(sum(t>mint&t<maxt),1)...
                                jj*j*ones(sum(t>mint&t<maxt),1) ii*ones(sum(t>mint&t<maxt),1)...
                                j*ones(sum(t>mint&t<maxt),1) i*ones(sum(t>mint&t<maxt),1)]];
                            pass_numspks = [pass_numspks; [length(spks) p jj*j ii j i]];    
                            field_vars = [field_vars; [field_var jj*j ii j i]];
                            pass_vel = [pass_vel;[vel(t>mint&t<maxt,1) p*ones(sum(t>mint&t<maxt),1) jj*j*ones(sum(t>mint&t<maxt),1)...
                                ii*ones(sum(t>mint&t<maxt),1) j*ones(sum(t>mint&t<maxt),1) i*ones(sum(t>mint&t<maxt),1)]];  
                            spk_vel = [spk_vel; [vel(t>min(field_ts(field_ts>mint))&t<max(field_ts(field_ts<maxt)),1)...
                                p*ones(sum(t>min(field_ts(field_ts>mint))&t<max(field_ts(field_ts<maxt))),1)...
                                jj*j*ones(sum(t>min(field_ts(field_ts>mint))&t<max(field_ts(field_ts<maxt))),1)...
                                ii*ones(sum(t>min(field_ts(field_ts>mint))&t<max(field_ts(field_ts<maxt))),1)...
                                j*ones(sum(t>min(field_ts(field_ts>mint))&t<max(field_ts(field_ts<maxt))),1)...
                                i*ones(sum(t>min(field_ts(field_ts>mint))&t<max(field_ts(field_ts<maxt))),1)]];
                            
                            if retro == 1
                                cmodes = [cmodes;[1 p jj*j ii j i]];
                                rImX = [rImX ImX(:,tt>mint&tt<maxt)];
                                rImX_phase = [rImX_phase;passImX_phase(tt>mint&tt<maxt)];
                                rspkImX = [rspkImX ImX(:,tt>min(field_ts(field_ts>mint))&tt<max(field_ts(field_ts<maxt)))];
                                rV = [rV;v(tt>min(field_ts(field_ts>mint))&tt<max(field_ts(field_ts<maxt)))];
                            elseif pro == 1
                                cmodes = [cmodes;[-1 p jj*j ii j i]];
                                pImX = [pImX ImX(:,tt>mint&tt<maxt)];
                                pImX_phase = [pImX_phase;passImX_phase(tt>mint&tt<maxt)];
                                pspkImX = [pspkImX ImX(:,tt>min(field_ts(field_ts>mint))&tt<max(field_ts(field_ts<maxt)))];
                                pV = [pV;v(tt>min(field_ts(field_ts>mint))&tt<max(field_ts(field_ts<maxt)))];
                            else
                                cmodes = [cmodes;[nan p jj*j ii j i]];
                                aImX = [aImX ImX(:,tt>mint&tt<maxt)];
                                aImX_phase = [aImX_phase;passImX_phase(tt>mint&tt<maxt)];
                                aspkImX = [aspkImX ImX(:,tt>min(field_ts(field_ts>mint))&tt<max(field_ts(field_ts<maxt)))];
                                aV = [aV;v(tt>min(field_ts(field_ts>mint))&tt<max(field_ts(field_ts<maxt)))];
                            end
                            
                        elseif p==length(pass_times) && sum(field_ts>pass_times(p)) > 1 && abs(path(end)) > 2 && abs(path(find(t>pass_times(p),1,'first'))) > 2
                            
                            %get pass spikes
                            spks = field_spks(field_ts>pass_times(p));
                            
                            %determine coding mode
                            retro = 0; pro = 0;
                            if sum(spks>0)/length(spks) > 2/3
                                retro = 1;
                            elseif sum(spks<0)/length(spks) > 2/3
                                pro = 1;
                            end
                            
                            %get first spike, apply third rotation, and scale so that spikes fall between 0 and 1
                            fspk = min(spks);
                            R_mu = [cos(-fspk),-sin(-fspk);sin(-fspk),cos(-fspk)];
                            spks = (R_mu*[cos(spks),sin(spks)]')';
                            spks = atan2(spks(:,2),spks(:,1));
                            passImX_phase = (R_mu*[cos(ImX_phase),sin(ImX_phase)]')';
                            passImX_phase = atan2(passImX_phase(:,2),passImX_phase(:,1));
                            passpath = (R_mu*[cos(path),sin(path)]')';
                            passpath = atan2(passpath(:,2),passpath(:,1));
                            scale = 1/max(spks);
                            spks = scale*spks; passpath = scale*passpath;
                            passImX_phase = scale*passImX_phase;
                            
                            mint = t(find(t>pass_times(p),1,'first'));
                            maxt = t(find(t>pass_times(p),1,'last'));
                            
                            %process vars
                            pass_spks = [pass_spks; [spks p*ones(length(spks),1) jj*j*ones(length(spks),1) ii*ones(length(spks),1)...
                                j*ones(length(spks),1) i*ones(length(spks),1)]];
                            pass_ts = [pass_ts; [field_ts(field_ts>mint&field_ts<maxt) p*ones(sum(field_ts>mint&field_ts<maxt),1)...
                                jj*j*ones(sum(field_ts>mint&field_ts<maxt),1) ii*ones(sum(field_ts>mint&field_ts<maxt),1)...
                                j*ones(sum(field_ts>mint&field_ts<maxt),1) i*ones(sum(field_ts>mint&field_ts<maxt),1)]];
                            pass_ct = [pass_ct; [median(field_ts(field_ts>mint&field_ts<maxt)) p jj*j ii j i]];
                            pass_path = [pass_path; [passpath(t>mint&t<maxt) p*ones(sum(t>mint&t<maxt),1)...
                                jj*j*ones(sum(t>mint&t<maxt),1) ii*ones(sum(t>mint&t<maxt),1)...
                                j*ones(sum(t>mint&t<maxt),1) i*ones(sum(t>mint&t<maxt),1)]];
                            pass_t = [pass_t; [t(t>mint&t<maxt) p*ones(sum(t>mint&t<maxt),1)...
                                jj*j*ones(sum(t>mint&t<maxt),1) ii*ones(sum(t>mint&t<maxt),1)...
                                j*ones(sum(t>mint&t<maxt),1) i*ones(sum(t>mint&t<maxt),1)]];
                            pass_numspks = [pass_numspks; [length(spks) p jj*j ii j i]];    
                            field_vars = [field_vars; [field_var jj*j ii j i]];
                            pass_vel = [pass_vel;[vel(t>mint&t<maxt,1) p*ones(sum(t>mint&t<maxt),1) jj*j*ones(sum(t>mint&t<maxt),1)...
                                ii*ones(sum(t>mint&t<maxt),1) j*ones(sum(t>mint&t<maxt),1) i*ones(sum(t>mint&t<maxt),1)]];  
                            spk_vel = [spk_vel; [vel(t>min(field_ts(field_ts>mint))&t<max(field_ts(field_ts<maxt)),1)...
                                p*ones(sum(t>min(field_ts(field_ts>mint))&t<max(field_ts(field_ts<maxt))),1)...
                                jj*j*ones(sum(t>min(field_ts(field_ts>mint))&t<max(field_ts(field_ts<maxt))),1)...
                                ii*ones(sum(t>min(field_ts(field_ts>mint))&t<max(field_ts(field_ts<maxt))),1)...
                                j*ones(sum(t>min(field_ts(field_ts>mint))&t<max(field_ts(field_ts<maxt))),1)...
                                i*ones(sum(t>min(field_ts(field_ts>mint))&t<max(field_ts(field_ts<maxt))),1)]];
                            
                            if retro == 1
                                cmodes = [cmodes;[1 p jj*j ii j i]];
                                rImX = [rImX ImX(:,tt>mint&tt<maxt)];
                                rImX_phase = [rImX_phase;passImX_phase(tt>mint&tt<maxt)];
                                rspkImX = [rspkImX ImX(:,tt>min(field_ts(field_ts>mint))&tt<max(field_ts(field_ts<maxt)))];
                                rV = [rV;v(tt>min(field_ts(field_ts>mint))&tt<max(field_ts(field_ts<maxt)))];
                            elseif pro == 1
                                cmodes = [cmodes;[-1 p jj*j ii j i]];
                                pImX = [pImX ImX(:,tt>mint&tt<maxt)];
                                pImX_phase = [pImX_phase;passImX_phase(tt>mint&tt<maxt)];
                                pspkImX = [pspkImX ImX(:,tt>min(field_ts(field_ts>mint))&tt<max(field_ts(field_ts<maxt)))];
                                pV = [pV;v(tt>min(field_ts(field_ts>mint))&tt<max(field_ts(field_ts<maxt)))];
                            else
                                cmodes = [cmodes;[nan p jj*j ii j i]];
                                aImX = [aImX ImX(:,tt>mint&tt<maxt)];
                                aImX_phase = [aImX_phase;passImX_phase(tt>mint&tt<maxt)];
                                aspkImX = [aspkImX ImX(:,tt>min(field_ts(field_ts>mint))&tt<max(field_ts(field_ts<maxt)))];
                                aV = [aV;v(tt>min(field_ts(field_ts>mint))&tt<max(field_ts(field_ts<maxt)))];
                            end
                            
                        end
                    end
                end
            end
            %determine times between coding events in the same session
            if size(pass_ct,1) >= session_start_idx
                tdiffs = repmat(pass_ct(session_start_idx:end,1)',length(pass_ct(session_start_idx:end,1)),1);
                tdiffs = tdiffs - tdiffs';tdiffs(tdiffs==0) = 999;
                tdiffs = reshape(abs(triu(tdiffs,1))',1,size(tdiffs,1)*size(tdiffs,2));
                tdiffs(tdiffs==0) = [];tdiffs(tdiffs==999) = 0;
                Tdiffs = [Tdiffs tdiffs];
                %             %determine corresponding coding changes (-1 diff,nan amb,1 same)
                cdiffs = repmat(cmodes(session_start_idx:end,1)',length(pass_ct(session_start_idx:end,1)),1);
                cdiffs = cdiffs.*cdiffs';
                cdiffs = reshape(triu(cdiffs,1)',1,size(cdiffs,1)*size(cdiffs,2));
                cdiffs(cdiffs==0) = [];
                Cdiffs = [Cdiffs cdiffs];    
            end
            session_start_idx = size(pass_ct,1)+1;            
            kappa = 50;
            cutoff = 3*sqrt(1-besseli(1,kappa)/besseli(0,kappa)); %three standard deviations
            phasebins = linspace(-2,3,numbins);
            if length(pImX_phase) > 1000
                z=1
                tic
                for nb = 1:numbins
                    [pPFR(:,nb,sidx), pweights(:,nb,sidx)] = P_estimator_uv(pImX,pImX_phase,phasebins(nb),kappa,cutoff);
                end
                toc
            end
            clear pImX pImX_phase
            if length(rImX_phase) > 1000
                z=2
                tic
                for nb = 1:numbins
                    [rPFR(:,nb,sidx), rweights(:,nb,sidx)] = P_estimator_uv(rImX,rImX_phase,phasebins(nb),kappa,cutoff);
                end
                toc
            end
            clear rImX rImX_phase
            %use for loops to avoid memory errors caused by large matrices 
            if length(aImX_phase) > 1000
                z=3
                tic
                for nb = 1:numbins
                    [aPFR(:,nb,sidx), aweights(:,nb,sidx)] = P_estimator_uv(aImX,aImX_phase,phasebins(nb),kappa,cutoff);
                end
                toc
            end
            clear aImX aImX_phase
            sidx = sidx+1;            
        end 
    end
    
    %process data co-ocurring with spikes
    vbins = linspace(vmin,20,numvbins);
    if length(pspkImX) > 1000
        z=4
        tic
        outsum   = nansum(pspkImX,2);
        outsumW  = nansum(abs(pspkImX),2);
        pspkWPLI(:,midx) = (outsum./outsumW).^2;
        pspkweights(:,midx) = size(pspkImX,2)*ones(length(freqVec),1);
        [pVFR(:,:,midx), pvweights(:,:,midx)] = VFR_WPLI(pspkImX,pV,vbins);
        toc
    end
    clear pspkImX pV
    if length(rspkImX) > 1000
        z=5
        tic
        outsum   = nansum(rspkImX,2);
        outsumW  = nansum(abs(rspkImX),2);
        rspkWPLI(:,midx) = (outsum./outsumW).^2;
        rspkweights(:,midx) = size(rspkImX,2)*ones(length(freqVec),1);
        [rVFR(:,:,midx), rvweights(:,:,midx)] = VFR_WPLI(rspkImX,rV,vbins);
        toc
    end
    clear rspkImX rV
    if length(aspkImX) > 1000
        z=6
        tic
        outsum   = nansum(aspkImX,2);
        outsumW  = nansum(abs(aspkImX),2);
        aspkWPLI(:,midx) = (outsum./outsumW).^2;
        aspkweights(:,midx) = size(aspkImX,2)*ones(length(freqVec),1);
        [aVFR(:,:,midx), avweights(:,:,midx)] = VFR_WPLI(aspkImX,aV,vbins);
        toc
    end
    clear aspkImX aV
    
    gaPFR  = cat(3,gaPFR,aPFR);
    gpPFR  = cat(3,gpPFR,pPFR);
    grPFR  = cat(3,grPFR,rPFR);
    gaweights = cat(3,gaweights,aweights);
    gpweights = cat(3,gpweights,pweights);
    grweights = cat(3,grweights,rweights);
    aPFR = sum(aPFR.*aweights,3)./sum(aweights,3); clear aweights
    pPFR = sum(pPFR.*pweights,3)./sum(pweights,3); clear pweights
    rPFR = sum(rPFR.*rweights,3)./sum(rweights,3); clear rweights
    
    %saves data to struct
    data.freqVec = freqVec;
    data.pass_spks = pass_spks;
    data.pass_numspks = pass_numspks;
    data.pass_vel = pass_vel;
    data.spk_vel = spk_vel;
    data.pass_ts = pass_ts;
    data.pass_path = pass_path;
    data.pass_t = pass_t;
    data.field_vars = field_vars;
    data.cmodes = cmodes;
    data.pass_ct = pass_ct;
    data.Cdiffs = Cdiffs;
    data.Tdiffs = Tdiffs;
    data.aPFR = aPFR;
    data.rPFR = rPFR;
    data.pPFR = pPFR;
    data.aVFR = aVFR(:,:,midx);
    data.rVFR = rVFR(:,:,midx);
    data.pVFR = pVFR(:,:,midx);
    data.rspkWPLI = rspkWPLI(:,midx);
    data.pspkWPLI = pspkWPLI(:,midx);
    data.aspkWPLI = aspkWPLI(:,midx);
    filename = sprintf('%s%s%s%s',Mice{i},strcat('Pro_Retro_normalized','\'),'data.mat');
    save(filename,'-struct','data');
    clear data
    midx = midx+1;
    
    %create group data
    gpass_spks = [gpass_spks;pass_spks];
    gpass_numspks = [gpass_numspks;pass_numspks];
    gpass_vel = [gpass_vel;pass_vel];
    gspk_vel = [gspk_vel;spk_vel];
    gpass_ts = [gpass_ts;pass_ts];
    gpass_path = [gpass_path;pass_path];
    gpass_t = [gpass_t;pass_t];
    gfield_vars = [gfield_vars;field_vars];
    gcmodes = [gcmodes;cmodes];
    gpass_ct = [gpass_ct;pass_ct];
    gCdiffs = [gCdiffs Cdiffs];
    gTdiffs = [gTdiffs Tdiffs];
    
end

gaPFR = nansum(gaPFR.*gaweights,3)./nansum(gaweights,3); clear gaweights
gpPFR = nansum(gpPFR.*gpweights,3)./nansum(gpweights,3); clear gpweights
grPFR = nansum(grPFR.*grweights,3)./nansum(grweights,3); clear grweights
gaspkWPLI = nansum(aspkWPLI.*aspkweights,2)./nansum(aspkweights,2); clear aspkweights
gpspkWPLI = nansum(pspkWPLI.*pspkweights,2)./nansum(pspkweights,2); clear pspkweights
grspkWPLI = nansum(rspkWPLI.*rspkweights,2)./nansum(rspkweights,2); clear rspkweights
gaVFR = nansum(aVFR.*avweights,3)./nansum(avweights,3); clear avweights
gpVFR = nansum(pVFR.*pvweights,3)./nansum(pvweights,3); clear pvweights
grVFR = nansum(rVFR.*rvweights,3)./nansum(rvweights,3); clear rvweights

%saves data to struct
data.freqVec = freqVec;
data.pass_spks = gpass_spks;
data.pass_numspks = gpass_numspks;
data.pass_vel = gpass_vel;
data.spk_vel = gspk_vel;
data.pass_ts = gpass_ts;
data.pass_path = gpass_path;
data.pass_t = gpass_t;
data.field_vars = gfield_vars;
data.cmodes = gcmodes;
data.pass_ct = gpass_ct;
data.Cdiffs = gCdiffs;
data.Tdiffs = gTdiffs;
data.aPFR = gaPFR;
data.rPFR = grPFR;
data.pPFR = gpPFR;
data.aVFR = gaVFR;
data.rVFR = grVFR;
data.pVFR = gpVFR;
data.rspkWPLI = grspkWPLI;
data.pspkWPLI = gpspkWPLI;
data.aspkWPLI = gaspkWPLI;
filename = sprintf('%s%s%s%s',dd,strcat('\Pro_Retro_normalized','\'),'data.mat');
save(filename,'-struct','data');

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
function centre = findCentre(NE,NW,SW,SE)

% The centre will be at the point of interception by the corner diagonals
a = (NE(2)-SW(2))/(NE(1)-SW(1)); % Slope for the NE-SW diagonal
b = (SE(2)-NW(2))/(SE(1)-NW(1)); % Slope for the SE-NW diagonal
c = SW(2);
d = NW(2);
x = (d-c+a*SW(1)-b*NW(1))/(a-b); % X-coord of centre
y = a*(x-SW(1))+c; % Y-coord of centre
centre = [x,y];

function [P, W] = P_estimator_uv(ImX,ImXp,pb,kappa,cutoff)
% delta = angle(exp(1j*ImXp).*conj(exp(1j*pb*ones(size(ImXp)))));
delta = angle(exp(1j*ImXp)*conj(exp(1j*pb)));
ImX(:,abs(delta)>cutoff) = [];
delta(abs(delta)>cutoff) = [];
outsum = ImX*vm_kernel(delta,kappa);
outsumW = abs(ImX)*vm_kernel(delta,kappa);
P = (outsum./outsumW).^2;
W = ones(size(ImX))*vm_kernel(delta,kappa);

%Von Mises Kernel
function r = vm_kernel(x,kappa)

C = 1/(2*pi*besseli(0,kappa));
r = C * exp(kappa*cos(x));

function [VFR,W] = VFR_WPLI(ImX,vel,vbins)

invh = 0.5;
cutoff = 3*sqrt(1/invh);
VFR = zeros(size(ImX,1),length(vbins));
W = zeros(size(ImX,1),length(vbins));

for vb=1:length(vbins)      
    [VFR(:,vb), W(:,vb)] = P_estimator(ImX,vel,vbins(vb),invh,cutoff);    
end

% Estimate the mean P for one position value
function [p, w] = P_estimator(P,vel,vbin,invh,cutoff)
% edge-corrected kernel density estimator
P(:,abs(vel-vbin)>cutoff) = [];
vel(abs(vel-vbin)>cutoff) = [];
num = P*gaussian_kernel((vel-vbin),invh);
den = abs(P)*gaussian_kernel((vel-vbin),invh);
p = (num./den).^2;
w = ones(size(P))*gaussian_kernel((vel-vbin),invh);
                                                    
% Gaussian kernel for the rate calculation
function r = gaussian_kernel(x,invh)
r =  invh/sqrt(2*pi)*exp(-0.5*x.*x*invh^2);
