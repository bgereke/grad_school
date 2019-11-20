function [data] = decodewindows(sessions,cellist,data)

%this version finds the error and average TFR for each time window moving
%successively across the session
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
scale = .275; % cm per pixel
binsize = 3; % in cm
dt = .04;
step = .01;
x = [];
xL = [];
v = [];
a = [];
vt = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load position, compute velocity, acceleration for all sessions
for s = 1:length(sessions)
    vfile = [sessions{s},'VT1.nvt'];
    [x_temp,y,vt_temp] = readVideoData(vfile,scale);
    [x_temp,y] = axesRotater(x_temp,y);  %modify for 2D
    x_temp = x_temp - min(x_temp);
    %estimate and smooth velocity (set max to 120 cm/s)
    v_temp = findVelLinear(x_temp,vt_temp);   
    v_temp(v_temp>=120) = 0.5*(v_temp(circshift((v_temp>=120),-3))+v_temp(circshift((v_temp>=120),3)));
    v_temp = smooth(v_temp,15);

    %estimate acceleration
    a_temp = zeros(size(v_temp));
    for i = 2:1:(length(v_temp)-1)
        a_temp(i) = (abs(v_temp(i+1))-abs(v_temp(i-1)))/(vt_temp(i+1)-vt_temp(i-1)); 
    end
    a_temp(length(v_temp)) = (abs(v_temp(length(v_temp)))-abs(v_temp(length(v_temp)-1)))...
        /(vt_temp(length(vt_temp))-vt_temp(length(vt_temp)-1));
    a_temp = smooth(a_temp,15);

    %[xL_temp] = linearizedirection(x_temp,vt_temp,x_temp,vt_temp,v_temp); %consider reflipping leftward runs
    xL_temp = x_temp;
    if size(x_temp,1)==1
        x = [x x_temp];xL = [xL xL_temp];
    else
        x = [x;x_temp];xL = [xL;xL_temp];
    end
    if size(v_temp,1) == 1
        v = [v v_temp];a = [a a_temp];
    else
        v = [v; v_temp];a = [a; a_temp];
    end
    if size(vt_temp,1) == 1
        vt = [vt vt_temp];
    else
        vt = [vt; vt_temp];
    end
        
end
'stillworking3'
%get "measured position" in terms of bin for the time window
data.x = nan(size(data.w,1),1); %mean position for window
data.vel = nan(size(data.w,1),1); %mean velocity for window
data.acc = nan(size(data.w,1),1); %mean accerlation for window
data.wt = nan(size(data.w,1),1); %mean time for window
tLength = max(xL) - min(xL);
numBins = ceil(tLength / binsize);
data.xbins(:,1) = binsize*(1:numBins)+min(xL)-0.5*binsize; %location bins

for ww = 1:size(data.w,1)  
    xspan = xL(vt>data.w(ww,1) & vt<=data.w(ww,2));
    vspan = v(vt>data.w(ww,1) & vt<=data.w(ww,2));
    aspan = a(vt>data.w(ww,1) & vt<=data.w(ww,2));
    tspan = vt(vt>data.w(ww,1) & vt<=data.w(ww,2));
    if isempty(xspan)
        [xspan, yspan, vspan, aspan] = GetSpikePos(tstart+dt/2,xL,y,vt,v,a);
    end
    x1 = mean(xspan);
    xdiff = (data.xbins - x1).^2;
    [~,bin_idx] = min(xdiff);
    data.x(ww) = data.xbins(bin_idx);
    data.vel(ww) = mean(vspan);
    data.acc(ww) = mean(aspan);
    data.wt(ww) = mean(tspan);
end
'stillworking4'        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get cell ratemaps
i = 0;
fid=fopen(cellist);
while 1
       tline = fgetl(fid);
       if ~ischar(tline) 
           break
       end
       i = i+1;
       
       %load spikes from all sessions and combine
       tempspikes = [];
       for s = 1:length(sessions)
            tfile = [sessions{s},tline];       
            tempspikes_temp = loadSpikes(tfile);
            if size(tempspikes_temp,1) == 1
                tempspikes = [tempspikes tempspikes_temp];
            else
                tempspikes = [tempspikes;tempspikes_temp];
            end
       end
       
       %compute place fields + gamma triggered fields for both directions
       [spikes{i,1},spkidx{i,1},spksfgz{i,1}] = ...
           get_spike_times_within_specific_eeg_windows(data.w(:,1),data.w(:,2),tempspikes,data.wt,data.wsgzmu,data.wfgzmu); %to theta filter spikes and give gamma labels
       %spikes{i,1} = tempspikes; %to NOT theta filter the spikes
       [tempspkx,spky,spkv] = GetSpikePos(spikes{i,1},data.x,y,data.wt,data.vel,data.acc);
       %tempspkxL = linearizedirection(tempspkx,spikes{i,1},x,vt,v);
       %tempspkxL = getSpikePos(spikes{i,1},xL,y,t);
       [temprate] = ratemap_decode(data.x,tempspkx,data.vel,spkv,data.wsgzmu,data.wfgzmu, spksfgz{i,1},data.w,data.xbins); %consider using kernel density estimator
       data.cellrates_pos (:,i) = temprate(:,1); 
       data.cellrates_neg (:,i) = temprate(:,2);
       data.cellrates_pos_sg(:,i) = temprate(:,3);
       data.cellrates_neg_sg(:,i) = temprate(:,4);
       data.cellrates_pos_fg(:,i) = temprate(:,5);
       data.cellrates_neg_fg(:,i) = temprate(:,6);
end
'stillworking5'       
%save('cr','data.cellrates'); return;
numcells = size(data.cellrates_pos,2);
numbins = size(data.cellrates_pos,1);

%define probability densities for each sliding window
data.Px = probX(xL,binsize); %prior, optionally incorporate state transition matrix
data.Px = ones(numbins,1)/numbins;  %removes prior
data.Prx = cell(size(data.w,1),1); %likelihood
data.Prxg = cell(size(data.w,1),1); %likelihood
data.Pxr = cell(size(data.w,1),1); %posterior
data.Pxrg = cell(size(data.w,1),1); %posterior
data.HPxr = cell(size(data.w,1),1); %posterior entropy
data.HPxrg = cell(size(data.w,1),1); %posterior entropy
data.xhat = cell(size(data.w,1),1); %decoded position
data.xhatg = cell(size(data.w,1),1); %decoded position
data.xerr = cell(size(data.w,1),1); %position prediction error
data.xerrg = cell(size(data.w,1),1); %position prediction error
data.wxerr = cell(size(data.w,1),1); %mean prediciton error for window
data.wxerrg = cell(size(data.w,1),1); %mean prediciton error for window
data.wHPxr = cell(size(data.w,1),1); %mean posterior entropy for window
data.wHPxrg = cell(size(data.w,1),1); %mean posterior entropy for window
data.numsp = nan(size(data.w,1),1); %number of spikes for window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%iterate through windows
for ww = 1:size(data.w,1)    
    tstart = data.w(ww,1);
    %numsteps = (data.w(ww,2)-tstart)/step;
    numsteps = 1;   %use for cycle based measure
    dt = data.w(ww,2)-tstart;   %use for cycle based measure
    minspikes = 2;
    mincells = 2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %iterate time steps within each window
    data.xhat{ww} = nan(numsteps,3);
    data.xhatg{ww} = nan(numsteps,3);
    data.xerr{ww} = nan(numsteps,3);
    data.xerrg{ww} = nan(numsteps,3);
    data.Pxr{ww} = nan(numbins,numsteps);
    data.Pxrg{ww} = nan(numbins,numsteps);
    data.HPxr{ww} = nan(numsteps,1);
    data.HPxrg{ww} = nan(numsteps,1);
    data.Prx{ww} = zeros(numbins,numsteps);
    data.Prxg{ww} = zeros(numbins,numsteps);
    
    for t = 1:numsteps
                
        %get number of spikes for each cell within the time window
        r = zeros(1,numcells);
        for i = 1:numcells
             r(1,i) = length(find(spikes{i,1}>tstart & spikes{i,1}<=tstart+dt)); 
        end     
        numactivecells = length(find(r>0));        
        
       %get probabilty the animals was in each bin
        if sum(r) >= minspikes && numactivecells>=mincells
            data.numsp(ww) = sum(r);
             if data.vel(ww) >= 0
                 for i=1:numbins
                     Prix = zeros(numcells,1);
                     for j = 1:numcells
                         Prix(j) = ((dt*data.cellrates_pos(i,j))^r(1,j))/factorial(r(1,j));
                         Prix(j) = Prix(j) * exp(-dt*data.cellrates_pos(i,j));
                     end                     
                     data.Prx{ww}(i,t) = prod(Prix);
                     data.Pxr{ww}(i,t) =  data.Px(i) * data.Prx{ww}(i,t);
                     clear Prix;
                     if data.wsgzmu(ww) >= 0.25 && data.wfgzmu(ww) <= 0.5
                         Prix = zeros(numcells,1);
                         for j = 1:numcells
                             Prix(j) = ((dt*data.cellrates_pos_sg(i,j))^r(1,j))/factorial(r(1,j));
                             Prix(j) = Prix(j) * exp(-dt*data.cellrates_pos_sg(i,j));
                         end
                         data.Prxg{ww}(i,t) = prod(Prix);
                         data.Pxrg{ww}(i,t) =  data.Px(i) * data.Prxg{ww}(i,t);
                         clear Prix;
                     elseif data.wsgzmu(ww) <= 0.5 && data.wfgzmu(ww) >= 0.25
                         Prix = zeros(numcells,1);
                         for j = 1:numcells
                             Prix(j) = ((dt*data.cellrates_pos_fg(i,j))^r(1,j))/factorial(r(1,j));
                             Prix(j) = Prix(j) * exp(-dt*data.cellrates_pos_fg(i,j));
                         end
                         data.Prxg{ww}(i,t) = prod(Prix);
                         data.Pxrg{ww}(i,t) =  data.Px(i) * data.Prxg{ww}(i,t);
                         clear Prix;
                     end
                 end
             else
                 for i=1:numbins
                     Prix = zeros(numcells,1);
                     for j = 1:numcells
                         Prix(j) = ((dt*data.cellrates_neg(i,j))^r(1,j))/factorial(r(1,j));
                         Prix(j) = Prix(j) * exp(-dt*data.cellrates_neg(i,j));
                     end
                     data.Prx{ww}(i,t) = prod(Prix);
                     data.Pxr{ww}(i,t) =  data.Px(i) * data.Prx{ww}(i,t);
                     clear Prix;
                     if data.wsgzmu(ww) >= 0.25 && data.wfgzmu(ww) <= 0.5
                         Prix = zeros(numcells,1);
                         for j = 1:numcells
                             Prix(j) = ((dt*data.cellrates_neg_sg(i,j))^r(1,j))/factorial(r(1,j));
                             Prix(j) = Prix(j) * exp(-dt*data.cellrates_neg_sg(i,j));
                         end
                         data.Prxg{ww}(i,t) = prod(Prix);
                         data.Pxrg{ww}(i,t) =  data.Px(i) * data.Prxg{ww}(i,t);
                         clear Prix;
                     elseif data.wsgzmu(ww) <= 0.5 && data.wfgzmu(ww) >= 0.25
                         Prix = zeros(numcells,1);
                         for j = 1:numcells
                             Prix(j) = ((dt*data.cellrates_neg_fg(i,j))^r(1,j))/factorial(r(1,j));
                             Prix(j) = Prix(j) * exp(-dt*data.cellrates_neg_fg(i,j));
                         end
                         data.Prxg{ww}(i,t) = prod(Prix);
                         data.Pxrg{ww}(i,t) =  data.Px(i) * data.Prxg{ww}(i,t);
                         clear Prix;
                     end
                 end
             end
                 
            %convert to probability distribution
            data.Pxr{ww}(:,t) = data.Pxr{ww}(:,t)/nansum(data.Pxr{ww}(:,t));
            data.Pxrg{ww}(:,t) = data.Pxrg{ww}(:,t)/nansum(data.Pxrg{ww}(:,t)); 
            
            %compute MAP estimate
            [~,maxind] = max(data.Pxr{ww}(:,t));
            if ~isempty(maxind)
                data.xhat{ww}(t,1) = data.xbins(maxind);
            end
            [~,maxind] = max(data.Pxrg{ww}(:,t));
            if ~isempty(maxind)
                data.xhatg{ww}(t,1) = data.xbins(maxind);
            end
            
            %compute Bayes' least squares estimate            
            xmean = 0;
            for i = 1:numbins
                xmean = nansum([xmean  i*data.Pxr{ww}(i,t)]);
            end
            data.xhat{ww}(t,2) = binsize*xmean/nansum(data.Pxr{ww}(:,t));
            xmean = 0;
            for i = 1:numbins
                xmean = nansum([xmean  i*data.Pxrg{ww}(i,t)]);
            end
            data.xhatg{ww}(t,2) = binsize*xmean/nansum(data.Pxrg{ww}(:,t));
            
            %compute posterior median estimate
%             [Pxr_sort, xbin_idx] = sort(data.Pxr{ww}(:,t));
%             Pxr_sum = cumsum(Pxr_sort);
%             sort_idx = find(Pxr_sum >= 0.5, 1) %errors on empty matrix
%             xbin_idx(sort_idx)
%             data.xhat{ww}(t,3) = xbin_idx(sort_idx);
                        
            %compute posterior entropy
            xH = 0;
            for i = 1:numbins
                xH = nansum([xH  data.Pxr{ww}(i,t)*log2(data.Pxr{ww}(i,t))]);
            end
            data.HPxr{ww}(t) = -xH;
            xH = 0;
            for i = 1:numbins
                xH = nansum([xH  data.Pxrg{ww}(i,t)*log2(data.Pxrg{ww}(i,t))]);
            end
            data.HPxrg{ww}(t) = -xH;
            
        end
               
        %compute prediciton errors
        data.xerr{ww}(t,:) = data.xhat{ww}(t,:) - data.x(ww);
        if data.vel(ww) < 0
            data.xerr{ww}(t,:) = -1*data.xerr{ww}(t,:);
        end
        data.xerrg{ww}(t,:) = data.xhatg{ww}(t,:) - data.x(ww);
        if data.vel(ww) < 0
            data.xerrg{ww}(t,:) = -1*data.xerrg{ww}(t,:);
        end

        tbin(t) = tstart;

        tstart = tstart + step;
    end
    
    %compute mean window prediciton error, posterior variance and entropy
    data.wxerr{ww} = nanmean(data.xerr{ww},1);
    data.wHPxr{ww} = nanmean(data.HPxr{ww});
    data.wxerrg{ww} = nanmean(data.xerrg{ww},1);
    data.wHPxrg{ww} = nanmean(data.HPxrg{ww});
    
end
    
