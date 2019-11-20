function [fastresults,slowresults,xbins,cellrates,fPxn,sPxn] = decodesession2(cellist,eegfile)

%this version uses slwo and fast gamma windows and sees what the encoding
%error was during those times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
scale = .275; % cm per pixel
binsize = 1.5; % in cm
dt = .08;%ie size of window around the gamma peaks
zval = 3;%must be above zscore to count as peak
minspikes = 1;%must be this many spikes in the time window to analyuze it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load stuff
[x,y,t] = readVideoData('VT1.nvt',scale);
[x,y] = axesRotater(x,y);
x = x - min(x);
v = findVelLinear(x);
[xL] = linearizedirection(x,t,x,t,smooth(v,15));
[eeg,eegTS] = loadEeg8(eegfile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the gamma peaks
%first cut off ends of track to get rid of the highpower slow gamma stuff
[eegds,eegTSds] = eegEndofTrackCutter(eeg,eegTS,x,y,t,10,10);
[slowpeaks,~] = gammaWindows2(eegds,25,45,2000,dt,zval);
[fastpeaks,~] = gammaWindows2(eegds,75,95,2000,dt,zval);
%get start and stop times for the windows around the peaks
stimes = peakstowin(slowpeaks,eegds,eegTSds,dt);
ftimes = peakstowin(fastpeaks,eegds,eegTSds,dt);
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
       cellrates (:,i) = temprate; %binsxCell
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%go through all the windows
[fastresults,fPxn] = bayesrecon(ftimes,cellrates,spikes,xL,t,v,binsize,minspikes);
[slowresults,sPxn] = bayesrecon(stimes,cellrates,spikes,xL,t,v,binsize,minspikes);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xfin,Pxn] = bayesrecon(times,rates,spikes,xL,t,v,binsize,minspikes)
%get p(x)
[Px,xbins] = probX(xL,binsize);

numcells = size(rates,2);
numbins = size(rates,1);
numwindows = size(times,1);
xfin = nan(numwindows,4);
Pxn = zeros(numbins,numwindows);

for k = 1:numwindows
    dt = abs(times(k,2) - times(k,1));
    %get number of spikes for each cell within the time window (obseved
    %activity)
    n = zeros(1,numcells);
    for i = 1:numcells
         n(1,i) = length(find(spikes{i,1}>times(k,1) & spikes{i,1}<=times(k,2))); 
    end     
    %get probabilty that the animal was in each bin, for time window k, given n
    if sum(n) >= minspikes
        for i=1:numbins   
            for j = 1:numcells
                Pnx(j) = (dt*rates(i,j)^n(1,j))/factorial(n(1,j));
                Pnx(j) = Pnx(j) * exp(-dt*rates(i,j));
            end
            Pxn(i,k) =  Px(i) * prod(Pnx);
            clear Pnx;
        end 
        [~,maxind] = max(Pxn(:,k));
        if ~isempty(maxind)
            xfin(k,1) = maxind; %the predicted bin
        else 
        end
    end
    %get "measured position" in terms of bin for the time window..
    xspan = xL(t>times(k,1) & t<=times(k,2));    
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
    %xfin(k,2) = mean(xL(t>times(k,1) & t<=times(k,2) )  ); %weighted average of position
    %get velocity during the time window
    xfin(k,4) = mean(v(t>times(k,1) & t<=times(k,2)));
end

%loop around the values that occurred on the end
for i =1:k
    if (xfin(i,1)<15 && xfin(i,2) > numbins - 15)
        xfin(i,1) = xfin(i,1) + numbins;
    elseif (xfin(i,2)<15 && xfin(i,1) > numbins - 15)
        xfin(i,2) = xfin(i,2) + numbins;
    end
end

%get the errors
xfin(:,3) = xfin(:,1) - xfin(:,2);