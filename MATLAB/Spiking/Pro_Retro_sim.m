%generate trajectory
clear
cd('C:\Data\Pro_Retro_data\')
load('data.mat')
[vpd,vels] = hist(pass_vel(:,1),50);
vpd = vpd/sum(vpd);
vcpd = cumsum(vpd);
segl = 1;
T = 1500;
dt = 1/30;
t = 0:dt:T;
pos = zeros(length(t),1);
pos(1) = 0;
ns = 1;
ri = rand(1);
idx = find(ri>=vcpd,1,'first');
v = vels(idx);

for i = 2:length(t)
    
    if t(i) > ns*segl
        ns = ns+1;
        ri = rand(1);
        idx = find(vcpd>=ri,1,'first');
        v = vels(idx);
    end
    pos(i) = pos(i-1) + dt*v;  
    
end

D = 100;
phase = mod(2*pos/D,2*pi)-pi;
phase = round(10*phase)/10;

%generate spikes
numcells = 50;
numspks = 10;%floor(mean(pass_numspks(:,1)));
spk_cen = round(10*(2*pi*rand(numcells,1)-pi))/10;
first_ts = -2;last_ts = -0.5;
spk_ts = linspace(first_ts,last_ts,numspks);
spk_phase = nan(numspks*50,numcells);
spk_t = nan(numspks*50,numcells);

for nc = 1:numcells
    idx = find(phase==spk_cen(nc));
    t_cen = t(idx);
    t_cen(diff(t_cen)<5) = [];    
    for ts = 1:numspks
        for p = 1:length(t_cen)
            [~,idx] = min((t-(t_cen(p)+spk_ts(ts))).^2);
            spk_t(ts*p,nc) = t(idx);
            spk_phase(ts*p,nc) = phase(idx);
        end
    end
end

% plot(t,phase,'k');hold on
% scatter(spk_t(:,1),spk_phase(:,1),'.r')
% plot([t(1) t(end)],[circ_mean(spk_phase(:,1)) circ_mean(spk_phase(:,1))],'r')

%determine coding modes/times
means = circ_mean(spk_phase);
numa = zeros(numcells,1);
nump = zeros(numcells,1);
numr = zeros(numcells,1);
modes = [];

for nc = 1:numcells
    idx = find(phase==spk_cen(nc));
    t_cen = t(idx);
    t_cen(diff(t_cen)<5) = [];    
    for p = 1:length(t_cen)
        spks_p = spk_phase((spk_t(:,nc)-t_cen(p))<=first_ts,nc);
        deltas = angle(exp(1j*spks_p).*conj(exp(1j*means(nc)*ones(size(spks_p)))));
        if sum(deltas>0)/length(deltas)>=2/3
            numr(nc) = numr(nc) + 1;
            modes = [modes;1 median(spk_t((spk_t(:,nc)-t_cen(p))<=first_ts,nc))];
        elseif sum(deltas<0)/length(deltas)>=2/3
            nump(nc) = nump(nc) + 1;
            modes = [modes;-1 median(spk_t((spk_t(:,nc)-t_cen(p))<=first_ts,nc))];
        else
            numa(nc) = numa(nc) + 1;
            modes = [modes;nan nan];
        end
    end
end

%determine proportions and plot cpd
propsr = numr./(numr+nump+numa);
propsp = nump./(numr+nump+numa);
propse = (numr + nump)./(numr+nump+numa);

[rpd,p] = hist(propsr,10);
rpd = rpd/sum(rpd);
rcpd = cumsum(rpd);
[ppd,p] = hist(propsp,10);
ppd = ppd/sum(ppd);
pcpd = cumsum(ppd);
[epd,p] = hist(propse,10);
epd = epd/sum(epd);
ecpd = cumsum(epd);
plot(p,rcpd,'r');hold on 
plot(p,pcpd,'b')
plot(p,ecpd,'k')
xlabel('Proportion of runs');xlim([0 1])
ylabel('Cumulative probability')

%determine times between coding events
clear Tdiffs Cdiffs
Tdiffs = [];Cdiffs = [];
tdiffs = repmat(modes(:,2)',length(modes(:,2)),1);
tdiffs = tdiffs - tdiffs';tdiffs(tdiffs==0) = 999;
tdiffs = reshape(abs(triu(tdiffs,1))',1,size(tdiffs,1)*size(tdiffs,2));
tdiffs(tdiffs==0) = [];tdiffs(tdiffs==999) = 0;
Tdiffs = [Tdiffs tdiffs];
%determine corresponding coding changes (-1 diff,nan amb,1 same)
cdiffs = repmat(modes(:,1)',length(modes(:,1)),1);
cdiffs = cdiffs.*cdiffs';
cdiffs = reshape(triu(cdiffs,1)',1,size(cdiffs,1)*size(cdiffs,2));
cdiffs(cdiffs==0) = [];
Cdiffs = [Cdiffs cdiffs];

same = Tdiffs(Cdiffs == 1);
notsame = Tdiffs(Cdiffs == -1);
[Ns,Xs] = hist(same,0:0.2:50);
[Nns,Xns] = hist(notsame,0:0.2:50);
figure;plot(Xs,Ns,'k');hold on; plot(Xns,Nns,'--k');xlim([0 5])
xlabel('time (sec)');ylabel('# coding pairs')
legend('same','different')