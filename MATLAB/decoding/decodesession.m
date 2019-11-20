function [cellrates,xbins,Px,n,spikes,Pxn,xfin] = decodesession (cellist,eegfile)

%this version finds the error and average TFR for each time window moving
%successively across the session
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
scale = .275; % cm per pixel
binsize = 2; % in cm
dt = .25;
step = .1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load stuff
[x,y,t] = readVideoData('VT1.nvt',scale);
[x,y] = axesRotater(x,y);
x = x - min(x);
v = findVelLinear(x);
[xL] = linearizedirection(x,t,x,t,smooth(v,15));
[eeg,eegTS] = loadEeg8(eegfile);
zTFRs = zscore(TFR_frequency_band(eeg,2000,5,25,45));
zTFRf = zscore(TFR_frequency_band(eeg,2000,5,70,100));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%more parameters
tstart = min(t);
numsteps = (max(t)-tstart)/step;
minspikes = 7;
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
       spikes{i,1} = loadSpikes2(tline);
       tempspkx = getSpikePos(spikes{i,1},x,y,t);
       tempspkxL = linearizedirection(tempspkx,spikes{i,1},x,t,v);
       [temprate, xbins] = ratemap_decode(xL,tempspkxL,binsize);
       cellrates (:,i) = temprate;
end

numcells = size(cellrates,2);
numbins = size(cellrates,1);

%get p(x) and p(n)
Px = probX(xL,binsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%go through all the time steps
xfin = nan(numsteps,6);
for k = 1:numsteps
    tstart = tstart + step;
    %get number of spikes for each cell within the time windoW
    n = zeros(1,numcells);
    for i = 1:numcells
         n(1,i) = length(find(spikes{i,1}>tstart & spikes{i,1}<=tstart+dt)); 
    end     
    %get probabilty the animals was in each bin
    if sum(n) >= minspikes
        for i=1:numbins   
            for j = 1:numcells
                Pnx(j) = (dt*cellrates(i,j)^n(1,j))/factorial(n(1,j));
                Pnx(j) = Pnx(j) * exp(-dt*cellrates(i,j));
            end
            Pxn(i,k) =  Px(i) * prod(Pnx);
            clear Pnx;
        end 
        [~,maxind] = max(Pxn(:,k));
        if ~isempty(maxind)
            xfin(k,1) = maxind; %the predicted bin
        end
    end
    %get "measured position" in terms of bin for the time window..
    xspan = xL(t>tstart & t<=tstart+dt);    
    x1 = xspan(1);
    x2 = xspan(end);
    xdiff = (xbins - x1).^2;
    [~,b1] = min(xdiff);
    xdiff = (xbins - x2).^2;
    [~,b2] = min(xdiff);

    if abs(b2 - b1) > numbins /2 
        b2 = b2 + numbins;
    end
    xfin(k,2) = mean([b1 b2]); %the measured bin
    %get velocity during the time window
    xfin(k,4) = mean(v(t>tstart & t<=tstart+dt));
    %get average tfr values for each time window
    xfin(k,5) = mean(zTFRs(eegTS>tstart & eegTS<+tstart+dt));
    xfin(k,6) = mean(zTFRf(eegTS>tstart & eegTS<+tstart+dt));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop around the values that occurred on the end
for i =1:numsteps
    if (xfin(i,1)<15 && xfin(i,2) > numbins - 15)
        xfin(i,1) = xfin(i,1) + numbins;
    elseif (xfin(i,2)<15 && xfin(i,1) > numbins - 15)
        xfin(i,2) = xfin(i,2) + numbins;
    end
end

%get the errors
xfin(:,3) = xfin(:,1) - xfin(:,2);







