function [Pxn,actualxbin,xbins,tbin] = decodewindow (cellist,window)

%this version finds the error and average TFR for each time window moving
%successively across the session
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
scale = .275; % cm per pixel
binsize = 3; % in cm
dt = .04;
step = .01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load stuff
[x,y,t] = readVideoData('VT1.nvt',scale);
[x,y] = axesRotater(x,y);
x = x - min(x);
v = findVelLinear(x);
[xL] = linearizedirection(x,t,x,t,smooth(v,15));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%more parameters
tstart = window(1,1);
numsteps = (window(1,2)-tstart)/step;
minspikes = 1;
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

%get p(x) (first iteration of inference)
Px = probX(xL,binsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%go through all the time steps
xfin = nan(numsteps,6);
Pxn = zeros(numbins,numsteps);
for k = 1:numsteps
    %get velocity during the time window
    xfin(k,4) = mean(v(t>tstart & t<=tstart+dt));
    %optional get Px|v (potential second itertation of inference)
    [~,Px] = probXV(Px,xL,v,xfin(k,4),binsize);
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
    
    tbin(k) = tstart;
    
    tstart = tstart + step;
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
actualxbin = xfin(:,2);
%get the errors
xfin(:,3) = xfin(:,1) - xfin(:,2);

[eeg1,eegTS1] = loadeeg8('CSC1.ncs');
bpt1 = fftbandpass(eeg1,2000,5,6,11,12);
bps1 = fftbandpass(eeg1,2000,23,25,55,57);
bpf1 = fftbandpass(eeg1,2000,63,65,95,97);
[eeg2,eegTS2] = loadeeg8('CSC2.ncs');
bpt2 = fftbandpass(eeg2,2000,5,6,11,12);
bps2 = fftbandpass(eeg2,2000,23,25,55,57);
bpf2 = fftbandpass(eeg2,2000,63,65,95,97);
[eeg3,eegTS3] = loadeeg8('CSC3.ncs');
bpt3 = fftbandpass(eeg3,2000,5,6,11,12);
bps3 = fftbandpass(eeg3,2000,23,25,55,57);
bpf3 = fftbandpass(eeg3,2000,63,65,95,97);
[eeg4,eegTS4] = loadeeg8('CSC4.ncs');
bpt4 = fftbandpass(eeg4,2000,5,6,11,12);
bps4 = fftbandpass(eeg4,2000,23,25,55,57);
bpf4 = fftbandpass(eeg4,2000,63,65,95,97);
[eeg5,eegTS5] = loadeeg8('CSC5.ncs');
bpt5 = fftbandpass(eeg5,2000,5,6,11,12);
bps5 = fftbandpass(eeg5,2000,23,25,55,57);
bpf5 = fftbandpass(eeg5,2000,63,65,95,97);

figure;
subplot(3,5,1);hmap(Pxn,tbin,xbins); hold on; plot(tbin,xbins(round(xfin(:,2))),'r');
subplot(3,5,2);hmap(Pxn,tbin,xbins); hold on; plot(tbin,xbins(round(xfin(:,2))),'r');
subplot(3,5,3);hmap(Pxn,tbin,xbins); hold on; plot(tbin,xbins(round(xfin(:,2))),'r');
subplot(3,5,4);hmap(Pxn,tbin,xbins); hold on; plot(tbin,xbins(round(xfin(:,2))),'r');
subplot(3,5,5);hmap(Pxn,tbin,xbins); hold on; plot(tbin,xbins(round(xfin(:,2))),'r');


subplot(3,5,6);plot(eegTS1(eegTS1>window(1,1) & eegTS1<window(1,2)),eeg1(eegTS1>window(1,1) & eegTS1<window(1,2)));
subplot(3,5,6); hold on; plot(eegTS1(eegTS1>window(1,1) & eegTS1<window(1,2)),bpt1(eegTS1>window(1,1) & eegTS1<window(1,2)),'k');xlim([tbin(1)-step/2 tbin(end)+step/2]);
subplot(3,5,6); title('Tetrode 1');

subplot(3,5,7);plot(eegTS2(eegTS2>window(1,1) & eegTS2<window(1,2)),eeg2(eegTS2>window(1,1) & eegTS2<window(1,2)));
subplot(3,5,7); hold on; plot(eegTS2(eegTS2>window(1,1) & eegTS2<window(1,2)),bpt2(eegTS2>window(1,1) & eegTS2<window(1,2)),'k');xlim([tbin(1)-step/2 tbin(end)+step/2]);
subplot(3,5,7); title('Tetrode 2');

subplot(3,5,8);plot(eegTS3(eegTS3>window(1,1) & eegTS3<window(1,2)),eeg3(eegTS3>window(1,1) & eegTS3<window(1,2)));
subplot(3,5,8); hold on; plot(eegTS3(eegTS3>window(1,1) & eegTS3<window(1,2)),bpt3(eegTS3>window(1,1) & eegTS3<window(1,2)),'k');xlim([tbin(1)-step/2 tbin(end)+step/2]);
subplot(3,5,8); title('Tetrode 3');

subplot(3,5,9);plot(eegTS4(eegTS4>window(1,1) & eegTS4<window(1,2)),eeg4(eegTS4>window(1,1) & eegTS4<window(1,2)));
subplot(3,5,9); hold on; plot(eegTS4(eegTS4>window(1,1) & eegTS4<window(1,2)),bpt4(eegTS4>window(1,1) & eegTS4<window(1,2)),'k');xlim([tbin(1)-step/2 tbin(end)+step/2]);
subplot(3,5,9); title('Tetrode 4');

subplot(3,5,10);plot(eegTS5(eegTS5>window(1,1) & eegTS5<window(1,2)),eeg5(eegTS5>window(1,1) & eegTS5<window(1,2)));
subplot(3,5,10); hold on; plot(eegTS5(eegTS5>window(1,1) & eegTS5<window(1,2)),bpt5(eegTS5>window(1,1) & eegTS5<window(1,2)),'k');xlim([tbin(1)-step/2 tbin(end)+step/2]);
subplot(3,5,10); title('Tetrode 5');


subplot(3,5,11); hold on; plot(eegTS1(eegTS1>window(1,1) & eegTS1<window(1,2)),bps1(eegTS1>window(1,1) & eegTS1<window(1,2)))
subplot(3,5,11); hold on; plot(eegTS1(eegTS1>window(1,1) & eegTS1<window(1,2)),bpf1(eegTS1>window(1,1) & eegTS1<window(1,2)),'r');xlim([tbin(1)-step/2 tbin(end)+step/2]);

subplot(3,5,12); hold on; plot(eegTS2(eegTS2>window(1,1) & eegTS2<window(1,2)),bps2(eegTS2>window(1,1) & eegTS2<window(1,2)))
subplot(3,5,12); hold on; plot(eegTS2(eegTS2>window(1,1) & eegTS2<window(1,2)),bpf2(eegTS2>window(1,1) & eegTS2<window(1,2)),'r');xlim([tbin(1)-step/2 tbin(end)+step/2]);

subplot(3,5,13); hold on; plot(eegTS3(eegTS3>window(1,1) & eegTS3<window(1,2)),bps3(eegTS3>window(1,1) & eegTS3<window(1,2)))
subplot(3,5,13); hold on; plot(eegTS3(eegTS3>window(1,1) & eegTS3<window(1,2)),bpf3(eegTS3>window(1,1) & eegTS3<window(1,2)),'r');xlim([tbin(1)-step/2 tbin(end)+step/2]);

subplot(3,5,14); hold on; plot(eegTS4(eegTS4>window(1,1) & eegTS4<window(1,2)),bps4(eegTS4>window(1,1) & eegTS4<window(1,2)))
subplot(3,5,14); hold on; plot(eegTS4(eegTS4>window(1,1) & eegTS4<window(1,2)),bpf4(eegTS4>window(1,1) & eegTS4<window(1,2)),'r');xlim([tbin(1)-step/2 tbin(end)+step/2]);

subplot(3,5,15); hold on; plot(eegTS5(eegTS5>window(1,1) & eegTS5<window(1,2)),bps5(eegTS5>window(1,1) & eegTS5<window(1,2)))
subplot(3,5,15); hold on; plot(eegTS5(eegTS5>window(1,1) & eegTS5<window(1,2)),bpf5(eegTS5>window(1,1) & eegTS5<window(1,2)),'r');xlim([tbin(1)-step/2 tbin(end)+step/2]);




