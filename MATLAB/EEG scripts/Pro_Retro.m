function Pro_Retro(mice,freqVec)

D = 96.5; %track diameter in cm

%create directory to store files
    dd = cd; %data directory
    dirInfo = dir(cd);
    found = 0;
    for ss=1:size(dirInfo,1)
        if dirInfo(ss).isdir
            if strcmp(dirInfo(ss).name,strcat('Pro_Retro_data','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(dd,strcat('\Pro_Retro_data','\'));
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

%declare global mouse vars
gfmoms = [];
gbounds = [];
gcom = [];
gnfields = [];
grew_dist = []; %data,cell,session,day,mouse
gpass_spks = []; %data,pass,cell,session,day,mouse
gpass_ts = []; %data,pass,cell,session,day,mouse
gpass_thphase = []; %data,pass,cell,session,day,mouse
gpass_path = []; %data,pass,cell,session,day,mouse
gpass_t = []; %data,pass,cell,session,day,mouse
gpass_numspks = []; %data,pass,cell,session,day,mouse
gpass_tenter = [];
gpass_texit = [];
grate_enter = [];
grate_exit = [];
gvel_enter = [];
gvel_exit = [];
gfield_vars = []; %data,cell,session,day,mouse
gcspkrats = []; %data,pass,cell,session,day,mouse
gcvelrats = []; %data,pass,cell,session,day,mouse
gctrats = []; %data,pass,cell,session,day,mouse
gcrrats = []; %data,pass,cell,session,day,mouse
gpass_ct = []; %data,pass,cell,session,day,mouse
gpass_vel = []; %data,pass,cell,session,day,mouse
gpass_rate = [];
cn = 0;     %cell number

for i = 1:numMice
    disp(sprintf('%s%s','Reading data for: ',Mice{i}));
    
    %create directory to store files
    dirInfo = dir(Mice{i});
    found = 0;
    for ss=1:size(dirInfo,1)
        if dirInfo(ss).isdir
            if strcmp(dirInfo(ss).name,strcat('Pro_Retro_data','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(Mice{i},strcat('Pro_Retro_data','\'));
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
    
    %declars vars being processed by animal  
    fmoms = [];
    bounds = [];
    com = [];
    nfields = [];
    rew_dist = []; %data,cell,session,day,mouse
    pass_spks = []; %data,pass,cell,session,day,mouse
    pass_ts = []; %data,pass,cell,session,day,mouse
    pass_thphase = []; %data,pass,cell,session,day,mouse
    pass_path = []; %data,pass,cell,session,day,mouse
    pass_t = []; %data,pass,cell,session,day,mouse
    pass_numspks = []; %data,pass,cell,session,day,mouse
    pass_tenter = [];
    pass_texit = [];
    pass_rate = [];
    rate_enter = [];
    rate_exit = [];
    vel_enter = [];
    vel_exit = [];
    field_vars = []; %data,cell,session,day,mouse
    cspkrats = []; %data,pass,cell,session,day,mouse
    crrats = []; %data,pass,cell,session,day,mouse
    cvelrats = []; %data,pass,cell,session,day,mouse
    ctrats = []; %data,pass,cell,session,day,mouse
    pass_ct = []; %data,pass,cell,session,day,mouse
    pass_vel = []; %data,pass,cell,session,day,mouse
    
    % Perform all analyses for each date
    for j=1:numDates
        %set directory
        disp(Dates{j});
        newpath = strcat(Mice{i},Dates{j});
        cd(newpath);
%         inFile = 'infile - Copy.txt';
        inFile = 'inFile.txt';
        
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
        
%         % read the file names of the references for the channels to be used for
%         % averaging
%         refid = fopen(refList,'r');
%         jj = 1;
%         while ~feof(refid)
%             str = fgetl(refid);
%             refs(jj) = {str};
%             jj = jj+1;
%         end
%         
%         % read the file names of the channels to be used for averaging
%         avgid = fopen(ch4avg,'r');
%         jj = 1;
%         while ~feof(avgid)
%             str = fgetl(avgid);
%             avgs(jj) = {str};
%             jj = jj+1;
%         end
%         numrefs = jj-1;
        
        Time = cell(numsessions,1);
        Pos = cell(numsessions,1);
        Vel = cell(numsessions,1);
        
        %loop through sessions
        for ii = 1:numsessions
            disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));

            % Get position data
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
            %some days might only have one recorded LED (via Dylan effect)
            if sum(isnan(frontX))/length(frontX) > 0.5 || sum(isnan(backX))/length(backX)>0.5
                % Set the field selection for reading the video files. 1 = Add parameter, 0 = skip
                % parameter
                fieldSelect(1) = 1; % Timestamps
                fieldSelect(2) = 1; % Extracted X
                fieldSelect(3) = 1; % Extracted Y
                fieldSelect(4) = 0; % Extracted Angel
                fieldSelect(5) = 0; % Targets
                fieldSelect(6) = 0; % Points
                % Do we return header 1 = Yes, 0 = No.
                extractHead = 0;
                % 5 different extraction modes, see help file for Nlx2MatVt
                extractMod = 1; % Extract all data
               [Time{ii}, frontX, frontY] = Nlx2MatVT(file,fieldSelect,extractHead,extractMod);
               Time{ii} = Time{ii}'/1000000;
               ind = find(frontX == 0);               
               frontX(ind) = NaN;
               frontY(ind) = NaN;
               backX = frontX; 
               backY = frontY;
            end
            
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
%             [Pos{ii},~] = localfit_ad(Time{ii},unwrap(Pos{ii}),2,0.04,10);
%             [~,Vel{ii}] = localfit_ad(Time{ii},Pos{ii},2,0.175,5);
            
            Pos{ii} = unwrap(Pos{ii});
            span = 6/(max(Time{ii})-min(Time{ii}));
            Pos{ii} = smooth(Pos{ii},span,'loess');
            Vel{ii} = abs(diff(Pos{ii}))./diff(Time{ii});
            Vel{ii} = [Vel{ii}(1);Vel{ii}];
            span = 6/(max(Time{ii})-min(Time{ii}));
            Vel{ii} = exp(smooth(log(Vel{ii}),span,'loess'));
            
            Time{ii}(isnan(Vel{ii})) = [];Pos{ii}(isnan(Vel{ii})) = [];
            Vel{ii}(isnan(Vel{ii})) = [];
            Pos{ii} = reshape(wraptopi(Pos{ii}),length(Pos{ii}),1);
            Vel{ii} = reshape(Vel{ii},length(Vel{ii}),1);

        end
        mt = zeros(numsessions,1);
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
            mt(ii) = min(Time{ii}); Time{ii} = Time{ii}-mt(ii); %tt = tt-mt(ii); 
            
            if ii ~= 1
               cn = cn - numcells; 
            end
            
            %loop through cells
            for jj=1:numcells
                
                cn = cn + 1;
                
                disp(sprintf('%s%s','Cell: ',num2str(jj)));                
                tfile = [sessions{ii},cells{jj}];                
                [TS] = loadSpikes(tfile);
                TS = TS-mt(ii);
                [spk_phase] = spikePos(TS,Pos{ii},Time{ii}); 
                
                %compute place field
                binsize = 1; %cm
                mapAxis = linspace(-pi,pi,pi*D/binsize);
                kappa = 100; % Smoothing factor when calculating the ratemap
                minr = 1; minp = 0.1; %min rate and proportion of peak rate
                
                V = [];P = [];S = [];T = [];SV = [];
                for kk = 1:numsessions
                    V = [V; Vel{kk}];
                    P = [P Pos{kk}'];
                    T = [T; Time{kk}-min(Time{kk})]; 
                    tfile = [sessions{kk},cells{jj}];                
                    [spkt] = loadSpikes(tfile);
                    if kk <= ii
                        spkt = spkt - mt(kk);
                    else
                        spkt = spkt-min(Time{kk});
                    end
                    [spkp] = spikePos(spkt,Pos{kk},Time{kk}-min(Time{kk}));
                    S = [S;spkp];         
                    [spkv] = spikePos(spkt,Vel{kk},Time{kk}-min(Time{kk}));
                    SV = [SV;spkv];
                end
                %remove spikes within some distance of reward location
                pbins = linspace(-pi,pi,360);
                vmap = zeros(length(pbins),1);
                for p=1:length(pbins)
                    vmap(p) = mean(V(abs(angle(exp(1j*P).*conj(exp(1j*pbins(p)*ones(size(P))))))<=1/360*pi));
                end
                [~,im] = min(vmap);
%                 spk_temp = S(abs(angle(exp(1j*S).*conj(exp(1j*pbins(im)*ones(size(S))))))>0.4 & SV>vmin);
                spk_temp = S(abs(angle(exp(1j*S).*conj(exp(1j*pbins(im)*ones(size(S))))))>0.4);
                %compute map on remaining spikes
                [map,~,cidx] = circle_map(spk_temp,P,T',kappa,mapAxis);
                
                %determine # place fields
                numfields = 0;
                active = map > 0.2*max(map) & map > 0.5;
                lcross = find(diff(active)==1)+1;
                rcross = find(diff(active)==-1)+1;
                if active(1) == 1 && active(end) == 0
                    lcross = [1 lcross];
                elseif active(1) == 0 && active(end) == 1
                    rcross = [rcross length(map)];
                end
                for l = 1:length(lcross)-1
                    rcross_l = min(rcross(rcross>lcross(l)));
                    if rcross_l - lcross(l) > 13
                        numfields = numfields + 1;
                    end
                end
                if max(rcross) > max(lcross)
                    if max(rcross) - max(lcross) > 13
                        numfields = numfields + 1;
                    end
                elseif length(map) - max(lcross) + min(rcross) > 13       
                       numfields = numfields + 1; 
                end
                
                %do field detection on biggest peak and set boundaries ~Kevin's method
                minb = nan; maxb = nan; pad = 0.1;
                for b = 1:length(mapAxis)
                    if cidx - b > 0
                        bb = cidx - b;
                    else
                        bb = length(mapAxis) - b + cidx;
                    end
                    if map(bb)<minr && map(bb)<minp*map(cidx)
                        minb = mapAxis(bb);
                        break;
                    end
                end
                for b = 1:length(mapAxis)
                    if cidx + b <= length(mapAxis)
                        bb = cidx + b;
                    else
                        bb = cidx + b - length(mapAxis);
                    end
                    if map(bb)<minr && map(bb)<minp*map(cidx)
                        maxb = mapAxis(bb);
                        break;
                    end
                end
                if minb > -pi+pad
                    minb = minb - pad;
                else
                    minb = minb - pad + 2*pi;
                end
                if maxb < pi-pad
                    maxb = maxb + pad;
                else
                    maxb = maxb + pad - 2*pi;
                end
                field_id = ones(length(spk_phase),1);
                fmap = map;
                if minb < maxb
                    field_id(spk_phase<minb) = 0;
                    field_id(spk_phase>maxb) = 0;
                    fmap(mapAxis<minb) = 0;
                    fmap(mapAxis>maxb) = 0;
                else
                    field_id(spk_phase<minb & spk_phase>maxb) = 0;
                    fmap(mapAxis<minb & mapAxis>maxb) = 0;
                end
                
                field_spks = spk_phase(field_id == 1);
                field_ts = TS(field_id == 1);
                nrminb = minb; %unrotated field boundaries
                nrmaxb = maxb; 
                fmean = circ_mean(mapAxis',fmap); %unrotated field mean
                f_mu = [cos(-fmean),-sin(-fmean);sin(-fmean),cos(-fmean)];
                raxis = (f_mu*[cos(mapAxis'),sin(mapAxis')]')';
                raxis = atan2(raxis(:,2),raxis(:,1));
                fstd = sqrt(sum(fmap.*raxis.^2)/sum(fmap)); %unrotated field width/standard deviation
                
                %plot and save field detection figs
                plot(Time{ii},Pos{ii},'-k'); hold on
                plot(field_ts,field_spks,'.r');title(strcat('Cell ',num2str(jj),'-',cells{jj}))
                plot([Time{ii}(1) Time{ii}(end)],[maxb maxb],'-b');
                plot([Time{ii}(1) Time{ii}(end)],[minb minb],'-b');
                if ~isempty(field_spks)
                    plot([Time{ii}(1) Time{ii}(end)],[circ_median(field_spks) circ_median(field_spks)],'-r');
                end
                plot([Time{ii}(1) Time{ii}(end)],[pbins(im) pbins(im)],'--b');
                xlabel('time (sec)');ylabel('track angle');   hold off             
                bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('PlaceFieldAnglesPR','\'),cells{jj}(1:end-2),'.bmp');
                %        epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('PlaceFieldAngles','\'),cells{jj}(1:end-2),'.eps');
                f = getframe(gcf);
                [pic, cmap] = frame2im(f);
                imwrite(pic,bmpImage,'bmp');
                
                %min distance of field to reward location
                rd = min([abs(angle(exp(1j*minb).*conj(exp(1j*pbins(im))))),abs(angle(exp(1j*maxb).*conj(exp(1j*pbins(im)))))]);
                rew_dist = [rew_dist;[rd jj ii j i]]; %data,cell,session,day,mouse
                
                %find field center
                field_mu = circ_mean(field_spks);                
                
                %rotate spikes to be centered at 0
                R_mu = [cos(-field_mu),-sin(-field_mu);sin(-field_mu),cos(-field_mu)];
                field_spks = (R_mu*[cos(field_spks),sin(field_spks)]')';
                field_spks = atan2(field_spks(:,2),field_spks(:,1));
                %                 plot(field_ts,field_spks,'.k')
                %                 ylim([-pi pi])
                %                 pause
                %rotate animal trajectory accordingly and break into single passes
                path = (R_mu*[cos(Pos{ii}');sin(Pos{ii}')])';
                path = atan2(path(:,2),path(:,1));   
                
                minb = (R_mu*[cos(minb);sin(minb)])';
                minb = atan2(minb(:,2),minb(:,1));
                maxb = (R_mu*[cos(maxb);sin(maxb)])';
                maxb = atan2(maxb(:,2),maxb(:,1));
                spk_phase = (R_mu*[cos(spk_phase),sin(spk_phase)]')';
                spk_phase = atan2(spk_phase(:,2),spk_phase(:,1));
                inspks = field_spks;
                ints = field_ts;
                
                %                 ImX_phase = (R_mu*[cos(ImX_phase),sin(ImX_phase)]')';
                %                 ImX_phase = atan2(ImX_phase(:,2),ImX_phase(:,1));
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
                vmin = 0.1; 
                num_passes = 0; %num passes that meet criteria
                pass_idx = zeros(size(pass_times)); %index of accepted passes             
                
                for p = 1:length(pass_times)+1
                    if p==1
                        numsp = sum(field_ts<pass_times(p));
                        if numsp>=1
                            mints = min(field_ts(field_ts<pass_times(p)));
                            maxts = max(field_ts(field_ts<pass_times(p)));
                            minv = min(Vel{ii}(Time{ii}>=mints&Time{ii}<=maxts));
                            meanv = mean(Vel{ii}(Time{ii}>=mints&Time{ii}<=maxts));
                            if meanv > vmin 
                                num_passes = num_passes + 1;
                                pass_idx(p) = p;
                            end
                        else
                            field_spks(field_ts<pass_times(p)) = nan;
                            field_ts(field_ts<pass_times(p)) = nan;
                        end
                    elseif p>1 && p<=length(pass_times)
                        numsp = sum(field_ts>pass_times(p-1)&field_ts<pass_times(p));
                        if numsp>=1
                            mints = min(field_ts(field_ts>pass_times(p-1)&field_ts<pass_times(p)));
                            maxts = max(field_ts(field_ts>pass_times(p-1)&field_ts<pass_times(p)));
                            minv = min(Vel{ii}(Time{ii}>=mints&Time{ii}<=maxts));
                            meanv = mean(Vel{ii}(Time{ii}>=mints&Time{ii}<=maxts));
                            if meanv > vmin 
                                num_passes = num_passes + 1;
                                pass_idx(p) = p;
                            end
                        else
                            field_spks(field_ts>pass_times(p-1)&field_ts<pass_times(p)) = nan;
                            field_ts(field_ts>pass_times(p-1)&field_ts<pass_times(p)) = nan;
                        end
                    elseif p>length(pass_times)    
                        numsp = sum(field_ts>pass_times(p-1));
                        if numsp>=1
                            mints = min(field_ts(field_ts>pass_times(p-1)));
                            maxts = max(field_ts(field_ts>pass_times(p-1)));
                            minv = min(Vel{ii}(Time{ii}>=mints&Time{ii}<=maxts));
                            meanv = mean(Vel{ii}(Time{ii}>=mints&Time{ii}<=maxts));
                            if meanv > vmin 
                                num_passes = num_passes + 1;
                                pass_idx(p) = p;
                            end
                        else
                            field_spks(field_ts>pass_times(p-1)) = nan;
                            field_ts(field_ts>pass_times(p-1)) = nan;
                        end
                    end
                end
                
                field_spks(isnan(field_spks)) = [];
                field_ts(isnan(field_ts)) = [];
                pass_idx(pass_idx == 0) = [];
                
                if num_passes > 3
                    %apply second rotation based on remaining spikes
%                     kappa = 150/(max(field_spks)-min(field_spks));
%                     mapdef = mapAxis(abs(mapAxis)<=pi/2);
%                     [passmap,~,~] = circle_map(field_spks,path',Time{ii}',kappa,mapdef);
%                     field_mu = circ_mean(mapdef',passmap);
                    field_mu = circ_median(field_spks);
%                     field_mu = circ_mean(first_spks);
%                     [field_var, ~] = circ_std(field_spks);
                    field_var = max(field_spks)-min(field_spks);
                    R_mu = [cos(-field_mu),-sin(-field_mu);sin(-field_mu),cos(-field_mu)];
                    field_spks = (R_mu*[cos(field_spks),sin(field_spks)]')';
                    field_spks = atan2(field_spks(:,2),field_spks(:,1));
                    %                     ImX_phase = (R_mu*[cos(ImX_phase),sin(ImX_phase)]')';
                    %                     ImX_phase = atan2(ImX_phase(:,2),ImX_phase(:,1));
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
%                     
%                     %plot single cell for SfN poster
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
                        if p==1 && sum(field_ts<pass_times(p)) > 1 && nanmean(cos(path(1:20))) < 0.5 ...
                                && sum(diff(Time{ii}(path>min(field_spks)&path<0&Time{ii}<pass_times(p)))) > 0 &&...
                                sum(diff(Time{ii}(path<max(field_spks)&path>0&Time{ii}<pass_times(p)))) > 0                      
                            
                            pass_spks = [pass_spks; [field_spks(field_ts<pass_times(p))...
                                p*ones(sum(field_ts<pass_times(p)),1) cn*ones(sum(field_ts<pass_times(p)),1)...
                                ii*ones(sum(field_ts<pass_times(p)),1) j*ones(sum(field_ts<pass_times(p)),1)...
                                i*ones(sum(field_ts<pass_times(p)),1)]];

                            pass_ts = [pass_ts; [field_ts(field_ts<pass_times(p)) p*ones(sum(field_ts<pass_times(p)),1)...
                                cn*ones(sum(field_ts<pass_times(p)),1) ii*ones(sum(field_ts<pass_times(p)),1)...
                                j*ones(sum(field_ts<pass_times(p)),1) i*ones(sum(field_ts<pass_times(p)),1)]];
                                                       
                            pass_ct = [pass_ct; [median(field_ts(field_ts<pass_times(p))) p cn ii j i rd]];
                            
                            pass_path = [pass_path; [path(Time{ii}<pass_times(p)) p*ones(sum(Time{ii}<pass_times(p)),1)...
                                cn*ones(sum(Time{ii}<pass_times(p)),1) ii*ones(sum(Time{ii}<pass_times(p)),1)...
                                j*ones(sum(Time{ii}<pass_times(p)),1) i*ones(sum(Time{ii}<pass_times(p)),1)]];

                            pass_t = [pass_t; [Time{ii}(Time{ii}<pass_times(p)) p*ones(sum(Time{ii}<pass_times(p)),1)...
                                cn*ones(sum(Time{ii}<pass_times(p)),1) ii*ones(sum(Time{ii}<pass_times(p)),1)...
                                j*ones(sum(Time{ii}<pass_times(p)),1) i*ones(sum(Time{ii}<pass_times(p)),1)]];
                            
                            pass_vel = [pass_vel; [Vel{ii}(Time{ii}<pass_times(p)) p*ones(sum(Time{ii}<pass_times(p)),1)...
                                cn*ones(sum(Time{ii}<pass_times(p)),1) ii*ones(sum(Time{ii}<pass_times(p)),1)...
                                j*ones(sum(Time{ii}<pass_times(p)),1) i*ones(sum(Time{ii}<pass_times(p)),1)]];

                            pass_numspks = [pass_numspks; [sum(field_ts<pass_times(p)) p cn ii j i]];
                            pass_tenter = [pass_tenter; [sum(diff(Time{ii}(path>min(field_spks)&path<0&Time{ii}<pass_times(p)))) p cn ii j i]];
                            pass_texit = [pass_texit; [sum(diff(Time{ii}(path<max(field_spks)&path>0&Time{ii}<pass_times(p)))) p cn ii j i]];
                            rate_enter = [rate_enter; [sum(field_spks(field_ts<pass_times(p))<0)/pass_tenter(end,1) p cn ii j i]];
                            rate_exit = [rate_exit; [sum(field_spks(field_ts<pass_times(p))>0)/pass_texit(end,1) p cn ii j i]];
                            vel_enter = [vel_enter; [abs(min(field_spks))/pass_tenter(end,1) p cn ii j i]];
                            vel_exit = [vel_exit; [max(field_spks)/pass_texit(end,1) p cn ii j i]];

                            cspkrats = [cspkrats;[sum(field_spks(field_ts<pass_times(p))>0)/pass_numspks(end,1) p cn ii j i]];
                            crrats = [crrats;[rate_exit(end,1)/(rate_enter(end,1)+rate_exit(end,1)) p cn ii j i]];
                            ctrats = [ctrats;[pass_texit(end,1)/(pass_tenter(end,1)+pass_texit(end,1)) p cn ii j i]];                           
                            
%                             cmidx = round(crrats(end,1)*256);
%                             if cmidx == 0
%                                 cmidx = 1;
%                             elseif cmidx > 256
%                                 cmidx = 256;
%                             end
%                             plot(field_ts(field_ts<pass_times(p))+tm-max(Time{ii}),field_spks(field_ts<pass_times(p)),'o','MarkerFaceColor',cm(cmidx,:),'MarkerEdgeColor','k');                                     
                            
                            %rate-based definition
                            %compute single pass spatial rate profile
                            mapdef = mapAxis(abs(mapAxis)<=pi/2);
                            [passmap,~,~] = circle_map(field_spks(field_ts<pass_times(p)),path(Time{ii}<pass_times(p))',Time{ii}(Time{ii}<pass_times(p))',kappa,mapdef);
                            com = [com; [circ_mean(mapdef',passmap) p cn ii j i]];
                            nfields = [nfields; [numfields p cn ii j i]];
                            pass_rate = [pass_rate;passmap'];
                            bounds = [bounds;nrminb nrmaxb];
                            fmoms = [fmoms;fmean fstd];
                            
                            if vel_exit(end,1)/(vel_enter(end,1)+vel_exit(end,1))>2/3
                                %                                 vcmodes = [vcmodes;[1 p cn ii j i]];
                                cvelrats = [cvelrats;[vel_exit(end,1)/(vel_enter(end,1)+vel_exit(end,1)) p cn ii j i]];
                            elseif vel_enter(end,1)/(vel_enter(end,1)+vel_exit(end,1))>2/3
                                %                                 vcmodes = [vcmodes;[-1 p cn ii j i]];
                                cvelrats = [cvelrats;[vel_exit(end,1)/(vel_enter(end,1)+vel_exit(end,1)) p cn ii j i]];
                            else
                                %                                 vcmodes = [vcmodes;[nan p cn ii j i]];
                                cvelrats = [cvelrats;[vel_exit(end,1)/(vel_enter(end,1)+vel_exit(end,1)) p cn ii j i]];
                            end                            
                            
                            mint = Time{ii}(Time{ii}<pass_times(p));mint = mint(1);
                            maxt = Time{ii}(Time{ii}<pass_times(p));maxt = maxt(end);
                            
                            field_vars = [field_vars; [field_var cn ii j i]];

                        elseif p>1 && p<=length(pass_times) && sum(field_ts>pass_times(p-1)&field_ts<pass_times(p)) > 1 ...
                                && sum(diff(Time{ii}(path>min(field_spks)&path<0&Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)))) > 0 ...
                            && sum(diff(Time{ii}(path<max(field_spks)&path>0&Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)))) > 0
                            
                            lap = p;
                            if nanmean(cos(path(1:20))) >= 0.5
                                lap = p-1;
                            end
                            
                            pass_spks = [pass_spks; [field_spks(field_ts>pass_times(p-1)&field_ts<pass_times(p))...
                                lap*ones(sum(field_ts>pass_times(p-1)&field_ts<pass_times(p)),1) cn*ones(sum(field_ts>pass_times(p-1)&field_ts<pass_times(p)),1)...
                                ii*ones(sum(field_ts>pass_times(p-1)&field_ts<pass_times(p)),1) j*ones(sum(field_ts>pass_times(p-1)&field_ts<pass_times(p)),1)...
                                i*ones(sum(field_ts>pass_times(p-1)&field_ts<pass_times(p)),1)]];

                            pass_ts = [pass_ts; [field_ts(field_ts>pass_times(p-1)&field_ts<pass_times(p)) lap*ones(sum(field_ts>pass_times(p-1)&field_ts<pass_times(p)),1)...
                                cn*ones(sum(field_ts>pass_times(p-1)&field_ts<pass_times(p)),1) ii*ones(sum(field_ts>pass_times(p-1)&field_ts<pass_times(p)),1)...
                                j*ones(sum(field_ts>pass_times(p-1)&field_ts<pass_times(p)),1) i*ones(sum(field_ts>pass_times(p-1)&field_ts<pass_times(p)),1)]];
                           
                            pass_ct = [pass_ct; [median(field_ts(field_ts>pass_times(p-1)&field_ts<pass_times(p))) p cn ii j i rd]];
                            pass_path = [pass_path; [path(Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)) lap*ones(sum(Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)),1)...
                                cn*ones(sum(Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)),1) ii*ones(sum(Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)),1)...
                                j*ones(sum(Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)),1) i*ones(sum(Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)),1)]];

                            pass_t = [pass_t; [Time{ii}(Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)) lap*ones(sum(Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)),1)...
                                cn*ones(sum(Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)),1) ii*ones(sum(Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)),1)...
                                j*ones(sum(Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)),1) i*ones(sum(Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)),1)]];
                            pass_vel = [pass_vel; [Vel{ii}(Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)) lap*ones(sum(Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)),1)...
                                cn*ones(sum(Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)),1) ii*ones(sum(Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)),1)...
                                j*ones(sum(Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)),1) i*ones(sum(Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)),1)]];

                            pass_numspks = [pass_numspks; [sum(field_ts>pass_times(p-1)&field_ts<pass_times(p)) lap cn ii j i]];  
                            pass_tenter = [pass_tenter; [sum(diff(Time{ii}(path>min(field_spks)&path<0&Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)))) lap cn ii j i]];
                            pass_texit = [pass_texit; [sum(diff(Time{ii}(path<max(field_spks)&path>0&Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p)))) lap cn ii j i]];
                            rate_enter = [rate_enter; [sum(field_spks(field_ts>pass_times(p-1)&field_ts<pass_times(p))<0)/pass_tenter(end,1) lap cn ii j i]];
                            rate_exit = [rate_exit; [sum(field_spks(field_ts>pass_times(p-1)&field_ts<pass_times(p))>0)/pass_texit(end,1) lap cn ii j i]];
                            vel_enter = [vel_enter; [abs(min(field_spks))/pass_tenter(end,1) lap cn ii j i]];
                            vel_exit = [vel_exit; [max(field_spks)/pass_texit(end,1) lap cn ii j i]];

                            cspkrats = [cspkrats;[sum(field_spks(field_ts>pass_times(p-1)&field_ts<pass_times(p))>0)/pass_numspks(end,1) lap cn ii j i]];
                            crrats = [crrats;[rate_exit(end,1)/(rate_enter(end,1)+rate_exit(end,1)) lap cn ii j i]];
                            ctrats = [ctrats;[pass_texit(end,1)/(pass_tenter(end,1)+pass_texit(end,1)) lap cn ii j i]];
                            
%                             cmidx = round(crrats(end,1)*256);
%                             if cmidx == 0
%                                 cmidx = 1;
%                             elseif cmidx > 256
%                                 cmidx = 256;
%                             end
%                             plot(field_ts(field_ts>pass_times(p-1)&field_ts<pass_times(p))+tm-max(Time{ii}),field_spks(field_ts>pass_times(p-1)&field_ts<pass_times(p)),'o','MarkerFaceColor',cm(cmidx,:),'MarkerEdgeColor','k');                           
                            
                            %rate-based definition
                            %compute single pass spatial rate profile
                            mapdef = mapAxis(abs(mapAxis)<=pi/2);
                            [passmap,~,~] = circle_map(field_spks(field_ts>pass_times(p-1)&field_ts<pass_times(p)),path(Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p))',...
                                Time{ii}(Time{ii}>pass_times(p-1)&Time{ii}<pass_times(p))',kappa,mapdef);
                            com = [com; [circ_mean(mapdef',passmap) p cn ii j i]];
                            nfields = [nfields; [numfields p cn ii j i]];
                            pass_rate = [pass_rate;passmap'];
                            
                            bounds = [bounds;nrminb nrmaxb];
                            fmoms = [fmoms;fmean fstd];

                            if vel_exit(end,1)/(vel_enter(end,1)+vel_exit(end,1))>2/3
%                                 vcmodes = [vcmodes;[1 lap cn ii j i]];
                                cvelrats = [cvelrats;[vel_exit(end,1)/(vel_enter(end,1)+vel_exit(end,1)) p cn ii j i]];
                            elseif vel_enter(end,1)/(vel_enter(end,1)+vel_exit(end,1))>2/3      
%                                 vcmodes = [vcmodes;[-1 lap cn ii j i]];
                                cvelrats = [cvelrats;[vel_exit(end,1)/(vel_enter(end,1)+vel_exit(end,1)) p cn ii j i]];
                            else
%                                 vcmodes = [vcmodes;[nan lap cn ii j i]];
                                cvelrats = [cvelrats;[vel_exit(end,1)/(vel_enter(end,1)+vel_exit(end,1)) p cn ii j i]];
                            end

                            mint = Time{ii}(Time{ii}>=pass_times(p-1)); mint = mint(1);
                            maxt = Time{ii}(Time{ii}<=pass_times(p)); maxt = maxt(end);
                            
                            field_vars = [field_vars; [field_var cn ii j i]];

                        elseif p>length(pass_times) && sum(field_ts>pass_times(p-1)) > 1 && nanmean(cos(path(end-20:end))) < 0.5 ...
                                && sum(diff(Time{ii}(path>min(field_spks)&path<0&Time{ii}>pass_times(p-1)))) > 0 && sum(diff(Time{ii}(path<max(field_spks)&path>0&Time{ii}>pass_times(p-1)))) > 0
                            
                            lap = p;
                            if nanmean(cos(path(1:20))) >= 0.5
                                lap = p-1;
                            end
                            
                            pass_spks = [pass_spks; [field_spks(field_ts>pass_times(p-1))...
                                lap*ones(sum(field_ts>pass_times(p-1)),1) cn*ones(sum(field_ts>pass_times(p-1)),1)...
                                ii*ones(sum(field_ts>pass_times(p-1)),1) j*ones(sum(field_ts>pass_times(p-1)),1)...
                                i*ones(sum(field_ts>pass_times(p-1)),1)]];

                            pass_ts = [pass_ts; [field_ts(field_ts>pass_times(p-1)) lap*ones(sum(field_ts>pass_times(p-1)),1)...
                                cn*ones(sum(field_ts>pass_times(p-1)),1) ii*ones(sum(field_ts>pass_times(p-1)),1)...
                                j*ones(sum(field_ts>pass_times(p-1)),1) i*ones(sum(field_ts>pass_times(p-1)),1)]];
                            
                            pass_ct = [pass_ct; [median(field_ts(field_ts>pass_times(p-1))) p cn ii j i rd]];
                            pass_path = [pass_path; [path(Time{ii}>pass_times(p-1)) lap*ones(sum(Time{ii}>pass_times(p-1)),1)...
                                cn*ones(sum(Time{ii}>pass_times(p-1)),1) ii*ones(sum(Time{ii}>pass_times(p-1)),1)...
                                j*ones(sum(Time{ii}>pass_times(p-1)),1) i*ones(sum(Time{ii}>pass_times(p-1)),1)]];

                            pass_t = [pass_t; [Time{ii}(Time{ii}>pass_times(p-1)) lap*ones(sum(Time{ii}>pass_times(p-1)),1)...
                                cn*ones(sum(Time{ii}>pass_times(p-1)),1) ii*ones(sum(Time{ii}>pass_times(p-1)),1)...
                                j*ones(sum(Time{ii}>pass_times(p-1)),1) i*ones(sum(Time{ii}>pass_times(p-1)),1)]];
                            pass_vel = [pass_vel; [Vel{ii}(Time{ii}>pass_times(p-1)) lap*ones(sum(Time{ii}>pass_times(p-1)),1)...
                                cn*ones(sum(Time{ii}>pass_times(p-1)),1) ii*ones(sum(Time{ii}>pass_times(p-1)),1)...
                                j*ones(sum(Time{ii}>pass_times(p-1)),1) i*ones(sum(Time{ii}>pass_times(p-1)),1)]];

                            pass_numspks = [pass_numspks; [sum(field_ts>pass_times(p-1)) p cn ii j i]];
                            pass_tenter = [pass_tenter; [sum(diff(Time{ii}(path>min(field_spks)&path<0&Time{ii}>pass_times(p-1)))) lap cn ii j i]];
                            pass_texit = [pass_texit; [sum(diff(Time{ii}(path<max(field_spks)&path>0&Time{ii}>pass_times(p-1)))) lap cn ii j i]];
                            rate_enter = [rate_enter; [sum(field_spks(field_ts>pass_times(p-1))<0)/pass_tenter(end,1) lap cn ii j i]];
                            rate_exit = [rate_exit; [sum(field_spks(field_ts>pass_times(p-1))>0)/pass_texit(end,1) lap cn ii j i]];
                            vel_enter = [vel_enter; [abs(min(field_spks))/pass_tenter(end,1) lap cn ii j i]];
                            vel_exit = [vel_exit; [max(field_spks)/pass_texit(end,1) lap cn ii j i]];

                            cspkrats = [cspkrats;[sum(field_spks(field_ts>pass_times(p-1))>0)/pass_numspks(end,1) lap cn ii j i]];
                            crrats = [crrats;[rate_exit(end,1)/(rate_enter(end,1)+rate_exit(end,1)) lap cn ii j i]];
                            ctrats = [ctrats;[pass_texit(end,1)/(pass_tenter(end,1)+pass_texit(end,1)) lap cn ii j i]];
  
%                             cmidx = round(crrats(end,1)*256);
%                             if cmidx == 0
%                                 cmidx = 1;
%                             elseif cmidx > 256
%                                 cmidx = 256;
%                             end
%                             plot(field_ts(field_ts>pass_times(p-1))+tm-max(Time{ii}),field_spks(field_ts>pass_times(p-1)),'o','MarkerFaceColor',cm(cmidx,:),'MarkerEdgeColor','k');
                          
                            %rate-based definition
                            %compute single pass spatial rate profile
                            mapdef = mapAxis(abs(mapAxis)<=pi/2);
                            [passmap,~,~] = circle_map(field_spks(field_ts>pass_times(p-1)),path(Time{ii}>pass_times(p-1))',Time{ii}(Time{ii}>pass_times(p-1))',kappa,mapdef);
                            com = [com; [circ_mean(mapdef',passmap) p cn ii j i]];
                            nfields = [nfields; [numfields p cn ii j i]];
                            pass_rate = [pass_rate;passmap'];
                            
                            bounds = [bounds;nrminb nrmaxb];
                            fmoms = [fmoms;fmean fstd];

                            if vel_exit(end,1)/(vel_enter(end,1)+vel_exit(end,1))>2/3
%                                 vcmodes = [vcmodes;[1 p cn ii j i]];
                                cvelrats = [cvelrats;[vel_exit(end,1)/(vel_enter(end,1)+vel_exit(end,1)) lap cn ii j i]];
                            elseif vel_enter(end,1)/(vel_enter(end,1)+vel_exit(end,1))>2/3      
%                                 vcmodes = [vcmodes;[-1 p cn ii j i]];
                                cvelrats = [cvelrats;[vel_exit(end,1)/(vel_enter(end,1)+vel_exit(end,1)) lap cn ii j i]];
                            else
%                                 vcmodes = [vcmodes;[nan p cn ii j i]];
                                cvelrats = [cvelrats;[vel_exit(end,1)/(vel_enter(end,1)+vel_exit(end,1)) lap cn ii j i]];
                            end
                           
                            mint = Time{ii}(Time{ii}>pass_times(p-1));mint = mint(1);
                            maxt = Time{ii}(Time{ii}>pass_times(p-1));maxt = maxt(end);
                            
                            field_vars = [field_vars; [field_var cn ii j i]];

                        end
                    end                    
                end                   
            end                
        end
    end
    
    %saves data to struct
    data.fmoms = fmoms;
    data.bounds = bounds;
    data.com = com;
    data.nfields = nfields;
    data.rew_dist = rew_dist;
    data.freqVec = freqVec;
    data.pass_spks = pass_spks; 
    data.pass_numspks = pass_numspks;
    data.pass_tenter = pass_tenter;
    data.pass_texit = pass_texit;
    data.rate_enter = rate_enter;
    data.rate_exit = rate_exit;
    data.vel_enter = vel_enter;
    data.vel_exit = vel_exit;
    data.pass_vel = pass_vel;
    data.pass_ts = pass_ts;
    data.pass_thphase = pass_thphase;
    data.pass_path = pass_path;
    data.pass_t = pass_t;
    data.pass_rate = pass_rate;
    data.field_vars = field_vars;
    data.cspkrats = cspkrats;
    data.crrats = crrats;
    data.cvelrats = cvelrats;
    data.ctrats = ctrats;
    data.pass_ct = pass_ct;
    filename = sprintf('%s%s%s%s',Mice{i},strcat('Pro_Retro_data','\'),'data.mat');
    save(filename,'-struct','data');
    clear data
    
    %create group data
    gfmoms = [gfmoms;fmoms];
    gbounds = [gbounds;bounds];
    gcom = [gcom;com];
    gnfields = [gnfields;nfields];
    grew_dist = [grew_dist;rew_dist];
    gpass_spks = [gpass_spks;pass_spks]; 
    gpass_numspks = [gpass_numspks;pass_numspks];
    gpass_tenter = [gpass_tenter;pass_tenter];
    gpass_texit = [gpass_texit;pass_texit];
    grate_enter = [grate_enter;rate_enter];
    grate_exit = [grate_exit;rate_exit];
    gvel_enter = [gvel_enter;vel_enter];
    gvel_exit = [gvel_exit;vel_exit];
    gpass_vel = [gpass_vel;pass_vel];
    gpass_ts = [gpass_ts;pass_ts];
    gpass_thphase = [gpass_thphase;pass_thphase];
    gpass_path = [gpass_path;pass_path];
    gpass_t = [gpass_t;pass_t];
    gpass_rate = [gpass_rate;pass_rate];
    gfield_vars = [gfield_vars;field_vars];
    gcspkrats = [gcspkrats;cspkrats];
    gcrrats = [gcrrats;crrats];
    gctrats = [gctrats;ctrats];
    gcvelrats = [gcvelrats;cvelrats];
    gpass_ct = [gpass_ct;pass_ct];
end


%saves data to struct
data.fmoms = gfmoms;
data.bounds = gbounds;
data.com = gcom;
data.nfields = gnfields;
data.rew_dist = grew_dist;
data.pass_spks = gpass_spks;
data.pass_numspks = gpass_numspks;
data.pass_tenter = gpass_tenter;
data.pass_texit = gpass_texit;
data.rate_enter = grate_enter;
data.rate_exit = grate_exit;
data.vel_enter = gvel_enter;
data.vel_exit = gvel_exit;
data.pass_vel = gpass_vel;
data.pass_ts = gpass_ts;
data.pass_thphase = gpass_thphase;
data.pass_path = gpass_path;
data.pass_t = gpass_t;
data.pass_rate = gpass_rate;
data.field_vars = gfield_vars;
data.cspkrats = gcspkrats;
data.crrats = gcrrats;
data.ctrats = gctrats;
data.cvelrats = gcvelrats;
data.pass_ct = gpass_ct;

filename = sprintf('%s%s%s%s',dd,strcat('\Pro_Retro_data','\'),'data.mat');
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

function [P, W] = P_estimator(ImX,ImXp,pb,kappa)
% edge-corrected kernel density estimator
delta = angle(repmat(exp(1j*ImXp),1,size(pb,2)).*conj(exp(1j*pb)));
outsum = ImX*vm_kernel(delta,kappa);
outsumW = abs(ImX)*vm_kernel(delta,kappa);
P = (outsum./outsumW).^2;
W = ones(size(ImX))*vm_kernel(delta,kappa);

function [P, W] = P_estimator_uv(ImX,ImXp,pb,kappa,cutoff)
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
    
    if frontX(ii)==0 || backX(ii)==0 || isnan(frontX(ii)) || isnan(backX(ii))
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
