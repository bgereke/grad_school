%% SECTION TITLE
% DESCRIPTIVE TEXT
%plot passes
[fvs,idx1] = sort(field_vars(:,1),'descend');
figure %all
n=0;
for i=1:length(fvs)
    n=n+1;
    varRows = ismember(pass_spks(:,2:6),cmodes(idx1(i),2:6),'rows');
    plot(pass_spks(varRows,1),n*ones(length(pass_spks(varRows,1)),1),'.k','MarkerSize',1);
    hold on;
end
hold on; plot([0 0],[0 n],'-r'); %xlim([-pi pi]);n
%% SECTION TITLE
% DESCRIPTIVE TEXT
figure %retro
n=0;
for i=1:length(fvs)
    if cmodes(idx1(i),1) == 1
        n=n+1;
        varRows = ismember(pass_spks(:,2:6),cmodes(idx1(i),2:6),'rows');
        plot(pass_spks(varRows,1),n*ones(length(pass_spks(varRows,1)),1),'.k','MarkerSize',1);
        hold on;
    end
end
hold on; plot([0 0],[0 n],'-r'); xlim([-pi pi]);n
%% SECTION TITLE
% DESCRIPTIVE TEXT
figure  %pro
n=0;
for i=1:length(fvs)
    if cmodes(idx1(i),1) == -1
        n=n+1;
        varRows = ismember(pass_spks(:,2:6),cmodes(idx1(i),2:6),'rows');
        plot(pass_spks(varRows,1),n*ones(length(pass_spks(varRows,1)),1),'.k','MarkerSize',1);
        hold on;
    end
end
hold on; plot([0 0],[0 n],'-r'); xlim([-pi pi]);n
%% SECTION TITLE
% DESCRIPTIVE TEXT
figure  %ambig
n=0;
for i=1:length(fvs)
    if isnan(cmodes(idx1(i),1))
        n=n+1;
        varRows = ismember(pass_spks(:,2:6),cmodes(idx1(i),2:6),'rows');
        plot(pass_spks(varRows,1),n*ones(length(pass_spks(varRows,1)),1),'.k','MarkerSize',1);
        hold on;
    end
end
hold on; plot([0 0],[0 n],'-r'); xlim([-pi pi]);n

%% SECTION TITLE
% DESCRIPTIVE TEXT
%cumulative probability distribution for coding events
cells = unique(cmodes(:,3));
props = zeros(size(cells));
props_r = zeros(size(cells));
props_p = zeros(size(cells));
props_sh = zeros(size(cells));
props_shr = zeros(size(cells));
props_shp = zeros(size(cells));
numshuf = 1;
for cid = 1:length(cells)
    pid = find(cmodes(:,3) == cells(cid)); %pass index
    props(cid) = sum(~isnan(cmodes(pid,1)))/length(pid);
    props_r(cid) = sum(cmodes(pid,1)==1)/length(pid);
    props_p(cid) = sum(cmodes(pid,1)==-1)/length(pid);
    num_retrosh= zeros(numshuf,length(pid));
    num_prosh = zeros(numshuf,length(pid));
    spk_bank = pass_spks(pass_spks(:,3)==cells(cid),1); %all the spikes from that cell
    pid = pid';
    for p = pid
        %shuffle spikes and recount
        for sh = 1:numshuf
            rndspks = randsample(spk_bank,pass_numspks(p,1));
            if sum(rndspks>0)/length(rndspks)>2/3
                num_retrosh(sh,p) = 1;
            elseif sum(rndspks<0)/length(rndspks)>2/3
                num_prosh(sh,p) = 1;
            end
        end
    end
    props_sh(cid) = (sum(sum(num_retrosh))+sum(sum(num_prosh)))/(numshuf*length(pid));    
    props_shr(cid) = sum(sum(num_retrosh))/(numshuf*length(pid));
    props_shp(cid) = sum(sum(num_prosh))/(numshuf*length(pid));
end

[pd,p] = hist(props,10);
[pd_r,p] = hist(props_r,10);
[pd_p,p] = hist(props_p,10);
pd = pd/sum(pd);
pd_r = pd_r/sum(pd_r);
pd_p = pd_p/sum(pd_p);
% props_sh = reshape(props_sh,1,size(props_sh,1)*size(props_sh,2));
% props_shr = reshape(props_shr,1,size(props_sh,1)*size(props_sh,2));
% props_shp = reshape(props_shp,1,size(props_sh,1)*size(props_sh,2));
[pd_sh,p] = hist(props_sh,10);
[pd_shr,p] = hist(props_shr,10);
[pd_shp,p] = hist(props_shp,10);
pd_sh = pd_sh/sum(pd_sh);
pd_shr = pd_shr/sum(pd_shr);
pd_shp = pd_shp/sum(pd_shp);
cpd = cumsum(pd);
cpd_r = cumsum(pd_r);
cpd_p = cumsum(pd_p);
cpd_sh = cumsum(pd_sh);
cpd_shr = cumsum(pd_shr);
cpd_shp = cumsum(pd_shp);
% figure;
% plot(p,pd,'k');hold on
% plot(p,pd_r,'r');
% plot(p,pd_p,'b');
% plot(p,pd_sh,'--k');
% plot(p,pd_shr,'--r');
% plot(p,pd_shp,'--b');
% legend('Data either','Data retro','Data pro','Shuffled either','Shuffled retro','Shuffled pro')
% xlabel('Proportion of runs');xlim([0 1])
% ylabel('Probability')
% title('Coding Mode Probability Distribution')
figure
plot(p,cpd,'k');hold on
plot(p,cpd_r,'r');
plot(p,cpd_p,'b');
plot(p,cpd_sh,'--k')
plot(p,cpd_shr,'--r')
plot(p,cpd_shp,'--b')
legend('Data either','Data retro','Data pro','Shuffled either','Shuffled retro','Shuffled pro')
xlabel('Proportion of runs');xlim([0 1])
ylabel('Cumulative probability')
title('Coding Mode Cumulative Probability Distribution')

%% SECTION TITLE
% DESCRIPTIVE TEXT
%hist version
same = Tdiffs(Cdiffs == 1);
notsame = Tdiffs(Cdiffs == -1);
[Ns,Xs] = hist(same,0:0.2:11);
[Nns,Xns] = hist(notsame,0:0.2:11);
figure;plot(Xs,Ns,'k');hold on; plot(Xns,Nns,'--k');xlim([0 5])
xlabel('time (sec)');ylabel('# coding pairs')
legend('same','different')
%hist version for velocity coding modes
same = Tdiffs(VCdiffs == 1);
notsame = Tdiffs(VCdiffs == -1);
[Ns,Xs] = hist(same,0:0.2:11);
[Nns,Xns] = hist(notsame,0:0.2:11);
figure;plot(Xs,Ns,'k');hold on; plot(Xns,Nns,'--k');xlim([0 5])
xlabel('time (sec)');ylabel('# coding pairs')
legend('same','different')
%% SECTION TITLE
% DESCRIPTIVE TEXT
%population cross correlation of ratios
td = 0:0.1:10;
sigma =0.5;
numshuf = 1000;
cv = zeros(length(td),numshuf);
cvsh = zeros(length(td),numshuf);
rat1 = tiedrank(Velratdiffs(:,1))-mean(tiedrank(Velratdiffs(:,1)));
rat2 = tiedrank(Velratdiffs(:,2))-mean(tiedrank(Velratdiffs(:,2)));
% rat1 = zRexitdiffs(:,1)-zRenterdiffs(:,1);
% rat2 = zRexitdiffs(:,2)-zRenterdiffs(:,2);
for i = 1:length(td)
    for s = 1:numshuf
        shidx1 = randsample(length(rat1),length(rat1),true);
        shidx2 = randsample(length(rat1),length(rat1),true);
        wave = exp(-0.5*(Tdiffs(shidx1)'-td(i)).^2/sigma^2);
        cv(i,s) = nansum(rat1(shidx1).*rat2(shidx1).*wave)/nansum(abs(rat1(shidx1).*rat2(shidx1)).*wave);
        cvsh(i,s) = nansum(rat1(shidx1).*rat2(shidx2).*wave)/nansum(abs(rat1(shidx1).*rat2(shidx2)).*wave);
        %    wave(abs(Tdiffs'-td(i) )>3*sigma) = [];
        %    rat1_temp = rat1; rat1_temp(abs(Tdiffs'-td(i))>3*sigma) = [];
        %    rat2_temp = rat2; rat2_temp(abs(Tdiffs'-td(i))>3*sigma) = [];
%         cv(i,s) = sum(rat1_temp(shidx).*rat2_temp(shidx).*wave)/sum(abs(rat1_temp(shidx).*rat2_temp(shidx)).*wave);
    end
end
figure;plot(td,median(cv,2),'r','LineWidth',1.5);hold on
plot(td,median(cvsh,2),'k','LineWidth',1.5);
legend('Real Data','Shuffled Data')
[cish] = prctile(cvsh,[2.5 97.5],2);
hold on; shadedplot(td, cish(:,2)', cish(:,1)',[0.7 0.7 0.7],[0.5 0.5 0.5]); alpha(0.8)
hold on;%xlim([0 5])
[ci] = prctile(cv,[2.5 97.5],2);
hold on; shadedplot(td, ci(:,2)', ci(:,1)',[1 0 0],[1 0.25 0.25]); alpha(0.5)
% hold on;plot([0 td(end)],[0 0],'--k')
xlabel('time between cell pairs (sec)');ylabel('Spearman Corr_R_I_,_R_I(time between cell pairs)')
box off

%% SECTION TITLE
% DESCRIPTIVE TEXT
%population partial cross correlation of ratios
td = 0:0.2:10;
sigma =0.5;
numshuf = 1000;
S = zeros(4); Ssh = zeros(4);
cv = zeros(length(td),numshuf);
cvsh = zeros(length(td),numshuf);
% rat5 = tiedrank(Rratdiffs(:,1))-mean(tiedrank(Rratdiffs(:,1)));
% rat6 = tiedrank(Rratdiffs(:,2))-mean(tiedrank(Rratdiffs(:,2)));
% rat3 = tiedrank(Spkratdiffs(:,1))-mean(tiedrank(Spkratdiffs(:,1)));
% rat4 = tiedrank(Spkratdiffs(:,2))-mean(tiedrank(Spkratdiffs(:,2)));
% rat1 = tiedrank(Velratdiffs(:,1))-mean(tiedrank(Velratdiffs(:,1)));
% rat2 = tiedrank(Velratdiffs(:,2))-mean(tiedrank(Velratdiffs(:,2)));
rat1 = tiedrank(Rratdiffs(:,1))-mean(tiedrank(Rratdiffs(:,1)));
rat2 = tiedrank(Rratdiffs(:,2))-mean(tiedrank(Rratdiffs(:,2)));
rat3 = tiedrank(Velratdiffs(:,1))-mean(tiedrank(Velratdiffs(:,1)));
rat4 = tiedrank(Velratdiffs(:,2))-mean(tiedrank(Velratdiffs(:,2)));
for i = 1:length(td)
    for s = 1:numshuf        
        shidx1 = randsample(length(rat1),length(rat1),true);
        shidx2 = randsample(length(rat1),length(rat1),true);
        wave = exp(-0.5*(Tdiffs(shidx1)'-td(i)).^2/sigma^2);
        S(1,1) = nansum(rat1(shidx1).*rat1(shidx1).*wave);
        S(2,2) = nansum(rat2(shidx1).*rat2(shidx1).*wave);
        S(3,3) = nansum(rat3(shidx1).*rat3(shidx1).*wave);
        S(4,4) = nansum(rat4(shidx1).*rat4(shidx1).*wave);
%         S(5,5) = sum(rat5(shidx1).*rat5(shidx1).*wave);
%         S(6,6) = sum(rat6(shidx1).*rat6(shidx1).*wave);
        S(1,2) = nansum(rat1(shidx1).*rat2(shidx1).*wave) ;
        S(1,3) = nansum(rat1(shidx1).*rat3(shidx1).*wave);
        S(1,4) = nansum(rat1(shidx1).*rat4(shidx1).*wave);
%         S(1,5) = sum(rat1(shidx1).*rat5(shidx1).*wave);
%         S(1,6) = sum(rat1(shidx1).*rat6(shidx1).*wave);
        S(2,3) = nansum(rat2(shidx1).*rat3(shidx1).*wave);
        S(2,4) = nansum(rat2(shidx1).*rat4(shidx1).*wave);
%         S(2,5) = sum(rat2(shidx1).*rat5(shidx1).*wave);
%         S(2,6) = sum(rat2(shidx1).*rat6(shidx1).*wave);
        S(3,4) = nansum(rat3(shidx1).*rat4(shidx1).*wave);
%         S(3,5) = sum(rat3(shidx1).*rat5(shidx1).*wave);
%         S(3,6) = sum(rat3(shidx1).*rat6(shidx1).*wave);
%         S(4,5) = sum(rat4(shidx1).*rat5(shidx1).*wave);
%         S(4,6) = sum(rat4(shidx1).*rat6(shidx1).*wave);
%         S(5,6) = sum(rat5(shidx1).*rat6(shidx1).*wave);
        S = triu(S) + triu(S,1)';
        Se = S(1:2,1:2)-S(1:2,3:4)*inv(S(3:4,3:4))*S(1:2,3:4)';
        R = diag(diag(Se).^-0.5)*Se*diag(diag(Se).^-0.5);
        cv(i,s) = R(2,1);
        %do same for shuffled
        Ssh(1,1) = nansum(rat1(shidx1).*rat1(shidx1).*wave);
        Ssh(2,2) = nansum(rat2(shidx2).*rat2(shidx2).*wave);
        Ssh(3,3) = nansum(rat3(shidx1).*rat3(shidx1).*wave);
        Ssh(4,4) = nansum(rat4(shidx2).*rat4(shidx2).*wave);
%         Ssh(5,5) = sum(rat5(shidx1).*rat5(shidx1).*wave);
%         Ssh(6,6) = sum(rat6(shidx2).*rat6(shidx2).*wave);
        Ssh(1,2) = nansum(rat1(shidx1).*rat2(shidx2).*wave);
        Ssh(1,3) = nansum(rat1(shidx1).*rat3(shidx1).*wave);
        Ssh(1,4) = nansum(rat1(shidx1).*rat4(shidx2).*wave);
%         Ssh(1,5) = sum(rat1(shidx1).*rat5(shidx1).*wave);
%         Ssh(1,6) = sum(rat1(shidx1).*rat6(shidx2).*wave);
        Ssh(2,3) = nansum(rat2(shidx2).*rat3(shidx1).*wave);
        Ssh(2,4) = nansum(rat2(shidx2).*rat4(shidx2).*wave);
%         Ssh(2,5) = sum(rat2(shidx2).*rat5(shidx1).*wave);
%         Ssh(2,6) = sum(rat2(shidx2).*rat6(shidx2).*wave);
        Ssh(3,4) = nansum(rat3(shidx1).*rat4(shidx2).*wave);
%         Ssh(3,5) = sum(rat3(shidx1).*rat5(shidx1).*wave);
%         Ssh(3,6) = sum(rat3(shidx1).*rat6(shidx2).*wave);
%         Ssh(4,5) = sum(rat4(shidx2).*rat5(shidx1).*wave);
%         Ssh(4,6) = sum(rat4(shidx2).*rat6(shidx2).*wave);
%         Ssh(5,6) = sum(rat5(shidx1).*rat6(shidx2).*wave);
        Ssh = triu(Ssh) + triu(Ssh,1)';
        Sesh = Ssh(1:2,1:2)-Ssh(1:2,3:4)*inv(Ssh(3:4,3:4))*Ssh(1:2,3:4)';
        Rsh = diag(diag(Sesh).^-0.5)*Sesh*diag(diag(Sesh).^-0.5);
        cvsh(i,s) = Rsh(2,1);
    end
    percent_done = i/length(td)*100
end
figure;plot(td,median(cv,2),'r','LineWidth',1.5);hold on
plot(td,median(cvsh,2),'k','LineWidth',1.5);
legend('Real Data','Shuffled Data')
[cish] = prctile(cvsh,[2.5 97.5],2);
hold on; shadedplot(td, cish(:,2)', cish(:,1)',[0.7 0.7 0.7],[0.5 0.5 0.5]); alpha(0.8)
hold on;%xlim([0 5])
[ci] = prctile(cv,[2.5 97.5],2);
hold on; shadedplot(td, ci(:,2)', ci(:,1)',[1 0 0],[1 0.25 0.25]); alpha(0.5)
% hold on;plot([0 td(end)],[0 0],'--k')
xlabel('time between cell pairs (sec)');ylabel('Spearman Corr_R_I_,_R_I(time between cell pairs)')
box off

% %time between coding pairs
% td = 0:0.02:10;
% invh = 5;
% same = Tdiffs(Cdiffs == 1);
% notsame = Tdiffs(Cdiffs == -1);
% Ns = zeros(size(td));
% Nns = zeros(size(td));
% for t = 1:length(td)
%     Ns(t) = sum(invh/sqrt(2*pi)*exp(-0.5*(same'-td(t)).*(same'-td(t))*invh^2));
%     Nns(t) = sum(invh/sqrt(2*pi)*exp(-0.5*(notsame'-td(t)).*(notsame'-td(t))*invh^2));
% end
% figure;plot(td,Ns/sum(Ns),'k');hold on; plot(td,Nns/sum(Nns),'--k');xlim([0 5])
% xlabel('Time between coding events from different cells (sec)')
% ylabel('Probability')
% 
% %same thing using matlab's ksdensity function
% same = Tdiffs(Cdiffs == 1);
% notsame = Tdiffs(Cdiffs == -1);
% [Ns,Xs] = ksdensity(same,0:0.2:11,'bandwidth',0.2);
% [Nns,Xns] = ksdensity(notsame,0:0.2:11,'bandwidth',0.2);
% figure;plot(Xs,Ns*length(same),'k');hold on; plot(Xns,Nns*length(notsame),'--k');xlim([0 5])

%% SECTION TITLE
% DESCRIPTIVE TEXT
%coding modes against lap # for session 1 w/bootstrapped confidence intervals
numlaps = 25;
nboot = 10000;
rprob = zeros(2,numlaps);
pprob = zeros(2,numlaps);
aprob = zeros(2,numlaps);
for lap = 1:numlaps
   
    numr = sum(cmodes(:,1)==1&cmodes(:,2)==lap&cmodes(:,4)==1);
    nump = sum(cmodes(:,1)==-1&cmodes(:,2)==lap&cmodes(:,4)==1);
    numa = sum(isnan(cmodes(:,1))&cmodes(:,2)==lap&cmodes(:,4)==1);
    
     if numr>1
        obsr = zeros(numr+nump+numa,1);       
        obsr(1:numr) = 1;        
        rprob(:,lap) = bootci(nboot,{@(x)[sum(x)/length(x)],obsr},'type','cper');        
     end
     if nump>1
         obsp = zeros(numr+nump+numa,1);
         obsp(1:nump) = 1;
         pprob(:,lap) = bootci(nboot,{@(x)[sum(x)/length(x)],obsp},'type','cper');
     end
     if numa>1
         obsa = zeros(numr+nump+numa,1);
         obsa(1:numa) = 1;   
         aprob(:,lap) = bootci(nboot,{@(x)[sum(x)/length(x)],obsa},'type','cper');
     end        
end
figure;
shadedplot(1:numlaps,rprob(1,:),rprob(2,:),[1 0.7 0.7], 'r'); hold on;
shadedplot(1:numlaps,pprob(1,:),pprob(2,:),[0.7 0.7 1],'b');hold on;
shadedplot(1:numlaps,aprob(1,:),aprob(2,:),[0.7 0.7 0.7],'k');hold on;
plot(1:numlaps,mean(rprob),'r','LineWidth',1.5);plot(1:numlaps,mean(pprob),'b','LineWidth',1.5);plot(1:numlaps,mean(aprob),'k','LineWidth',1.5)
xlabel('lap #');ylabel('Probability');
%legend('Retrospective','Prospective','Ambiguous')
figure;hist(cmodes(cmodes(:,4)==1,2),1:20);
xlabel('Lap number');ylabel('Count');

%% SECTION TITLE
% DESCRIPTIVE TEXT
%coding modes against time w/boostrapped confidence intervals
%have to run on data with t = t - min(t);
numbins = 25;
tbins = linspace(0,600,numbins);
nboot = 10000;     %number of bootstrap samples
rprob = zeros(2,numbins-1);
pprob = zeros(2,numbins-1);
aprob = zeros(2,numbins-1);
for t = 1:numbins-1
   
    numr = sum(cmodes(:,1)==1&pass_ct(:,1)>=tbins(t)&pass_ct(:,1)<tbins(t+1));
    nump = sum(cmodes(:,1)==-1&pass_ct(:,1)>=tbins(t)&pass_ct(:,1)<tbins(t+1));
    numa = sum(isnan(cmodes(:,1))&pass_ct(:,1)>=tbins(t)&pass_ct(:,1)<tbins(t+1));
    
    if numr>1
        obsr = zeros(numr+nump+numa,1);       
        obsr(1:numr) = 1;        
        rprob(:,t) = bootci(nboot,{@(x)[sum(x)/length(x)],obsr},'type','cper');        
     end
     if nump>1
         obsp = zeros(numr+nump+numa,1);
         obsp(1:nump) = 1;
         pprob(:,t) = bootci(nboot,{@(x)[sum(x)/length(x)],obsp},'type','cper');
     end
     if numa>1
         obsa = zeros(numr+nump+numa,1);
         obsa(1:numa) = 1;   
         aprob(:,t) = bootci(nboot,{@(x)[sum(x)/length(x)],obsa},'type','cper');
     end   
    
end

figure;
shadedplot(tbins(1:end-1),rprob(1,:),rprob(2,:),[1 0.7 0.7], 'r'); hold on;
shadedplot(tbins(1:end-1),pprob(1,:),pprob(2,:),[0.7 0.7 1],'b');hold on;
shadedplot(tbins(1:end-1),aprob(1,:),aprob(2,:),[0.7 0.7 0.7],'k');hold on;
plot(tbins(1:end-1),mean(rprob),'r','LineWidth',1.5);plot(tbins(1:end-1),mean(pprob),'b','LineWidth',1.5);plot(tbins(1:end-1),mean(aprob),'k','LineWidth',1.5)
xlabel('time (sec)');ylabel('Probability');
%legend('Retrospective','Prospective','Ambiguous')
figure;hist(pass_ct(:,1),numbins);
xlabel('times (sec)');ylabel('Count');

%% SECTION TITLE
% DESCRIPTIVE TEXT
%coding modes against running speed w/bootstrapped confidence intervals
cv = zeros(size(cmodes,1),1); %mean running speed during spiking
for c = 1:size(cmodes,1)
    idx = ismember(spk_vel(:,2:6),cmodes(c,2:6),'rows');
    cv(c) = nanmean(spk_vel(idx,1));   
end
dr=0.02;
rsp = 0:dr:0.6;         %runnning speed bins
nboot = 10000;     %number of bootstrap samples
rprob = zeros(2,length(rsp)-1);
pprob = zeros(2,length(rsp)-1);
aprob = zeros(2,length(rsp)-1);
for v = 1:length(rsp)-1
    
    numr = sum(cmodes(:,1)==1&cv>=rsp(v)&cv<rsp(v+1));
    nump = sum(cmodes(:,1)==-1&cv>=rsp(v)&cv<rsp(v+1));
    numa = sum(isnan(cmodes(:,1))&cv>=rsp(v)&cv<rsp(v+1)); 
    
    if numr>1
        obsr = zeros(numr+nump+numa,1);       
        obsr(1:numr) = 1;        
        rprob(:,v) = bootci(nboot,{@(x)[sum(x)/length(x)],obsr},'type','cper');        
     end
     if nump>1
         obsp = zeros(numr+nump+numa,1);
         obsp(1:nump) = 1;
         pprob(:,v) = bootci(nboot,{@(x)[sum(x)/length(x)],obsp},'type','cper');
     end
     if numa>1
         obsa = zeros(numr+nump+numa,1);
         obsa(1:numa) = 1;   
         aprob(:,v) = bootci(nboot,{@(x)[sum(x)/length(x)],obsa},'type','cper');
     end   
    
end

figure;
shadedplot(rsp(1:end-1)+dr/2,rprob(1,:),rprob(2,:),[1 0.7 0.7], 'r'); hold on;
shadedplot(rsp(1:end-1)+dr/2,pprob(1,:),pprob(2,:),[0.7 0.7 1],'b');hold on;
shadedplot(rsp(1:end-1)+dr/2,aprob(1,:),aprob(2,:),[0.7 0.7 0.7],'k');hold on;
plot(rsp(1:end-1)+dr/2,mean(rprob),'r','LineWidth',1.5);plot(rsp(1:end-1)+dr/2,mean(pprob),'b','LineWidth',1.5);plot(rsp(1:end-1)+dr/2,mean(aprob),'k','LineWidth',1.5)
xlabel('Running speed (cm/sec)');ylabel('Probability');
%legend('Retrospective','Prospective','Ambiguous')
figure;hist(cv,rsp);
xlabel('Running speed (cm/sec)');ylabel('Count');

%% SECTION TITLE
% DESCRIPTIVE TEXT
%Running Speed Against Coding Metrics

v = (vel_enter(:,1).*pass_tenter(:,1)+vel_exit(:,1).*pass_texit(:,1))./(pass_tenter(:,1)+pass_texit(:,1));
figure
plot(v,cspkrats(:,1),'.k')
xlabel('running speed')
ylabel('spike count proportions')
figure
plot(v,rate_exit(:,1)./(rate_exit(:,1)+rate_enter(:,1)),'.k')
xlabel('running speed')
ylabel('rate proportions')
figure
plot(v,zspk_exit(:,1)-zspk_enter(:,1),'.k')
xlabel('running speed')
ylabel('spike count z-score differences')
figure
plot(v,zrate_exit(:,1)-zrate_enter(:,1),'.k')
xlabel('running speed')
ylabel('rate z-score differences')

%% SECTION TITLE
% DESCRIPTIVE TEXT
%coding modes against running speed bootstrapped linear regression
cv = zeros(size(cmodes,1),1); %mean running speed during spiking
for c = 1:size(cmodes,1)
    idx = ismember(spk_vel(:,2:6),cmodes(c,2:6),'rows');
    cv(c) = mean(spk_vel(idx,1));    
end
dr=2;
rsp = 5:dr:25;         %runnning speed bins
nboot = 10000;     %number of bootstrap samples
rprobs = zeros(nboot,length(rsp)-1);
pprobs = zeros(nboot,length(rsp)-1);
aprobs = zeros(nboot,length(rsp)-1);

for v = 1:length(rsp)-1    
    numr = sum(cmodes(:,1)==1&cv>=rsp(v)&cv<rsp(v+1));
    nump = sum(cmodes(:,1)==-1&cv>=rsp(v)&cv<rsp(v+1));
    numa = sum(isnan(cmodes(:,1))&cv>=rsp(v)&cv<rsp(v+1)); 
    if numr>1
        obsr = zeros(numr+nump+numa,1);       
        obsr(1:numr) = 1;        
        rprobs(:,v) = bootstrp(nboot,@(x)[sum(x)/length(x)],obsr);        
     end
     if nump>1
         obsp = zeros(numr+nump+numa,1);
         obsp(1:nump) = 1;
         pprobs(:,v) = bootstrp(nboot,@(x)[sum(x)/length(x)],obsp);
     end
     if numa>1
         obsa = zeros(numr+nump+numa,1);
         obsa(1:numa) = 1;   
         aprobs(:,v) = bootstrp(nboot,@(x)[sum(x)/length(x)],obsa);
     end   
end

rstats = bootstrp(nboot,@myregress,repmat(rsp(1:end-1),nboot,1),rprobs);
pstats = bootstrp(nboot,@myregress,repmat(rsp(1:end-1),nboot,1),pprobs);
astats = bootstrp(nboot,@myregress,repmat(rsp(1:end-1),nboot,1),aprobs);

hist(rstats(:,1),100);
xlabel('slope');ylabel('count');title('retrospective');figure
hist(pstats(:,1),100);
xlabel('slope');ylabel('count');title('prospective');figure
hist(astats(:,1),100);
xlabel('slope');ylabel('count');title('ambiguous');

hist(rstats(:,2),50);figure
hist(pstats(:,2),50);figure
hist(astats(:,2),50);

hist(rstats(:,3),50);figure
hist(pstats(:,3),50);figure
hist(astats(:,3),50);

hist(rstats(:,4),50);figure
hist(pstats(:,4),50);figure
hist(astats(:,4),50);

hist(rstats(:,5),50);figure
hist(pstats(:,5),50);figure
hist(astats(:,5),50);

%% SECTION TITLE
% DESCRIPTIVE TEXT
%coding mode probabilities against pass spike number
dn = 2;
nspk = 5:dn:50;         %spike number bins
nboot = 10000;     %number of bootstrap samples
rprob = zeros(2,length(nspk)-1);
pprob = zeros(2,length(nspk)-1);
aprob = zeros(2,length(nspk)-1);
for n = 1:length(nspk)-1
    
    numr = sum(cmodes(:,1)==1&pass_numspks(:,1)>=nspk(n)&pass_numspks(:,1)<nspk(n+1));
    nump = sum(cmodes(:,1)==-1&pass_numspks(:,1)>=nspk(n)&pass_numspks(:,1)<nspk(n+1));
    numa = sum(isnan(cmodes(:,1))&pass_numspks(:,1)>=nspk(n)&pass_numspks(:,1)<nspk(n+1)); 
    
    if numr>1
        obsr = zeros(numr+nump+numa,1);       
        obsr(1:numr) = 1;        
        rprob(:,n) = bootci(nboot,{@(x)[sum(x)/length(x)],obsr},'type','cper');        
     end
     if nump>1
         obsp = zeros(numr+nump+numa,1);
         obsp(1:nump) = 1;
         pprob(:,n) = bootci(nboot,{@(x)[sum(x)/length(x)],obsp},'type','cper');
     end
     if numa>1
         obsa = zeros(numr+nump+numa,1);
         obsa(1:numa) = 1;   
         aprob(:,n) = bootci(nboot,{@(x)[sum(x)/length(x)],obsa},'type','cper');
     end   
    
end

figure;
shadedplot(nspk(1:end-1)+dn/2,rprob(1,:),rprob(2,:),[1 0.7 0.7], 'r'); hold on;
shadedplot(nspk(1:end-1)+dn/2,pprob(1,:),pprob(2,:),[0.7 0.7 1],'b');hold on;
shadedplot(nspk(1:end-1)+dn/2,aprob(1,:),aprob(2,:),[0.7 0.7 0.7],'k');hold on;
plot(nspk(1:end-1)+dn/2,mean(rprob),'r','LineWidth',1.5);plot(nspk(1:end-1)+dn/2,mean(pprob),'b','LineWidth',1.5);plot(nspk(1:end-1)+dn/2,mean(aprob),'k','LineWidth',1.5)
xlabel('# spikes');ylabel('Probability');
%legend('Retrospective','Prospective','Ambiguous')
figure;hist(pass_numspks(:,1),nspk);
xlabel('# spikes');ylabel('Count');

%% SECTION TITLE
% DESCRIPTIVE TEXT
%coding mode probabilities against field size
dn = 0.0025;
nspk = 0:dn:0.035;         %runnning speed bins
nboot = 10000;     %number of bootstrap samples
rprob = zeros(2,length(nspk)-1);
pprob = zeros(2,length(nspk)-1);
aprob = zeros(2,length(nspk)-1);
for n = 1:length(nspk)-1
    
    numr = sum(cmodes(:,1)==1&field_vars(:,1)>=nspk(n)&field_vars(:,1)<nspk(n+1));
    nump = sum(cmodes(:,1)==-1&field_vars(:,1)>=nspk(n)&field_vars(:,1)<nspk(n+1));
    numa = sum(isnan(cmodes(:,1))&field_vars(:,1)>=nspk(n)&field_vars(:,1)<nspk(n+1)); 
    
    if numr>1
        obsr = zeros(numr+nump+numa,1);       
        obsr(1:numr) = 1;        
        rprob(:,n) = bootci(nboot,{@(x)[sum(x)/length(x)],obsr},'type','cper');        
     end
     if nump>1
         obsp = zeros(numr+nump+numa,1);
         obsp(1:nump) = 1;
         pprob(:,n) = bootci(nboot,{@(x)[sum(x)/length(x)],obsp},'type','cper');
     end
     if numa>1
         obsa = zeros(numr+nump+numa,1);
         obsa(1:numa) = 1;   
         aprob(:,n) = bootci(nboot,{@(x)[sum(x)/length(x)],obsa},'type','cper');
     end   
    
end

figure;
shadedplot(nspk(1:end-1)+dn/2,rprob(1,:),rprob(2,:),[1 0.7 0.7], 'r'); hold on;
shadedplot(nspk(1:end-1)+dn/2,pprob(1,:),pprob(2,:),[0.7 0.7 1],'b');hold on;
shadedplot(nspk(1:end-1)+dn/2,aprob(1,:),aprob(2,:),[0.7 0.7 0.7],'k');hold on;
plot(nspk(1:end-1)+dn/2,mean(rprob),'r','LineWidth',1.5);plot(nspk(1:end-1)+dn/2,mean(pprob),'b','LineWidth',1.5);plot(nspk(1:end-1)+dn/2,mean(aprob),'k','LineWidth',1.5)
xlabel('field size');ylabel('Probability');
%legend('Retrospective','Prospective','Ambiguous')
figure;hist(field_vars(:,1),nspk);
xlabel('field size');ylabel('Count');


%% SECTION TITLE
% DESCRIPTIVE TEXT
%phase precession plots

%all
% figure
cells = unique(crrats(:,3));
slope = zeros(length(cells),1);
theta_0 = slope;
circ_corr = slope;
for nc = 1:length(cells)
   passes = find(crrats(:,3)==cells(nc));
   pb = repmat(linspace(-pi,pi,100),sum(pass_spks(:,3)==cells(nc)),1);
   delta = angle(repmat(exp(1j*pass_thphase(pass_spks(:,3)==cells(nc),1)),1,size(pb,2)).*conj(exp(1j*pb)));
   r = sum(exp(250*cos(delta)));   
   [~,idx] = max(r); 
   R_mu = [cos(-pb(1,idx)),-sin(-pb(1,idx));sin(-pb(1,idx)),cos(-pb(1,idx))];
   pos = [];phase = [];
   for p = 1:length(passes)
       rows = ismember(pass_spks(:,2:6),crrats(passes(p),2:6),'rows');
       pos = [pos;pass_spks(rows,1)];%-min(pass_spks(pass_spks(:,3)==cells(nc),1));
%        pos = pos/(max(pass_spks(pass_spks(:,3)==cells(nc),1))-min(pass_spks(pass_spks(:,3)==cells(nc),1)));
       phase_temp = (R_mu*[cos(pass_thphase(rows,1)),sin(pass_thphase(rows,1))]')';
       phase_temp = atan2(phase_temp(:,2),phase_temp(:,1));
       phase = [phase;phase_temp];
   end
   %        scatter([pos;pos],[phase;phase+2*pi],50,'.k')
   %        hold on
   [slope(nc), theta_0(nc), circ_corr(nc)] = circLinRegress(pos, phase);
%    pause
%    hold off
end
ppc = cells(circ_corr<-0.1);
% xlabel('normalized position in field');ylabel('theta phase');title('ambiguous')

%% SECTION TITLE
% DESCRIPTIVE TEXT
%single trial phase precession correlations, offsets, slopes, and ranges
cells = unique(crrats(:,3));
slope = [];
theta_0 = [];
circ_corr = [];
zspkdiffs = []; spkprops = [];
zrdiffs = []; rprops = [];
abins = -20:0.1:20;
R = zeros(size(abins));
for nc = 1:length(cells)
   passes = find(crrats(:,3)==cells(nc));
   pb = repmat(linspace(-pi,pi,100),sum(pass_spks(:,3)==cells(nc)),1);
   delta = angle(repmat(exp(1j*pass_thphase(pass_spks(:,3)==cells(nc),1)),1,size(pb,2)).*conj(exp(1j*pb)));
   r = sum(exp(250*cos(delta)));   
   [~,idx] = max(r); 
   R_mu = [cos(-pb(1,idx)),-sin(-pb(1,idx));sin(-pb(1,idx)),cos(-pb(1,idx))];
   Ph = [];Ps = [];
   for p = 1:length(passes)       
       rows = ismember(pass_spks(:,2:6),crrats(passes(p),2:6),'rows');
       spkprops = [spkprops;cspkrats(passes(p),1)];
       rprops = [rprops;rate_exit(passes(p))/(rate_exit(passes(p))+rate_enter(passes(p)))];
       zspkdiffs = [zspkdiffs;zspk_exit(passes(p),1) - zspk_enter(passes(p),1)];
       zrdiffs = [zrdiffs;zrate_exit(passes(p),1) - zrate_enter(passes(p),1)];
       pos = pass_spks(rows,1);
       phase = (R_mu*[cos(pass_thphase(rows,1)),sin(pass_thphase(rows,1))]')';
       phase = atan2(phase(:,2),phase(:,1));
%        t = pass_ts(rows,1);
%        if sum(diff(t)<0.025&diff(phase)>0)>0
%            tidx = diff(t)<0.025;
%            pidx = diff(phase)>0;
%            pos(find(tidx&pidx)+1) = [];
%            phase(find(tidx&pidx)+1) = [];
%            t(find(tidx&pidx)+1) = [];
%        end
       [slope_tmp, theta_0_tmp, circ_corr_tmp] = circLinRegress(pos, phase);      
       slope = [slope;slope_tmp];
       theta_0 = [theta_0;theta_0_tmp];
       circ_corr = [circ_corr;circ_corr_tmp];
%        scatter([pos;pos;pos],[phase;phase+2*pi;phase-2*pi],50,'.k')
%        figure(1)
%        scatter(pos,phase,50,'.k');hold on
%        scatter(pos,phase-2*pi,50,'.r')
%        scatter(pos,phase+2*pi,50,'.r');
%        title(['slope = ',num2str(slope(cnt)),' offset = ',num2str(theta_0(cnt)),' r = ',num2str(circ_corr(cnt))]);       
%        plot(pos,slope(cnt)*pos+theta_0(cnt),'-b')
%        hold off;
%        Ph = [Ph;phase];Ps = [Ps;pos]; 
%        figure(2)
%        for a = 1:length(abins)
%            R(a) = sqrt(nanmean(cos(phase-abins(a)*pos))^2+nanmean(sin(phase-abins(a)*pos))^2);
%        end
%        plot(abins,R);xlabel('slope');ylabel('R')
%        pause
   end
%    [slope, theta_0, circ_corr] = circLinRegress(Ps, Ph);
%    title(['slope = ',num2str(slope),' offset = ',num2str(theta_0),' r = ',num2str(circ_corr)]);
%    plot(Ps,slope*Ps+theta_0,'-b')   
%    hold off
%    
%    for a = 1:length(abins)
%        R(a) = sqrt(nanmean(cos(Ph-abins(a)*Ps))^2+nanmean(sin(Ph-abins(a)*Ps))^2);
%    end
%    figure(2)
%    plot(abins,R);xlabel('slope');ylabel('R')
%    pause
end
% numb = 10000;
% idx = reshape(randsample(length(rprops),length(rprops)*numb,true),length(rprops),numb);
% r = zeros(numb,1);p = zeros(numb,1);
% for i = 1:numb
%     [r(i),p(i)] = corr(rprops(idx(:,i)),circ_corr(idx(:,i)));
% end
% hist(r,100);xlabel('correletion of metric and theta phase-postion correlation');ylabel('# resamples')
% figure;hist(p,100);xlabel('p-values');ylabel('# resamples')
%% 

%shuffled ratio distributions
numshuf = 1000;
crates_she = zeros(length(cvelrats),numshuf);
crates_shx = zeros(length(cvelrats),numshuf);
crates_sh = zeros(length(cvelrats),numshuf);
crates = zeros(length(cvelrats),numshuf);
for cvid = 1:length(cvelrats)
    cidx = find(cvelrats(:,3)==cvelrats(cvid,3));
    spkbank = pass_spks(pass_spks(:,3)==cvelrats(cvid,3),1);
%     for sh = 1:numshuf
%         rndspks = randsample(spkbank,pass_numspks(cvid,1),true);
%         crates_she(cvid,sh) = sum(rndspks<0);
%         crates_shx(cvid,sh) = sum(rndspks>0);
%         crates_she(cvid,sh) = crates_she(cvid,sh)/randsample(pass_tenter(cidx),1);
%         crates_shx(cvid,sh) = crates_shx(cvid,sh)/randsample(pass_texit(cidx),1);       
% %         crates_she(cvid,sh) = sum(rndspks>0)/rndtexit/(sum(rndspks>0)/rndtexit+sum(rndspks<0)/rndtenter);
% %         crates_shx(cvid,sh) = sum(rndspks>0)/rndtexit/(sum(rndspks>0)/rndtexit+sum(rndspks<0)/rndtenter);
%     end
    rndenter = randsample(rate_enter(cidx,1).*pass_tenter(cidx,1),numshuf,true);
    rndexit = randsample(rate_exit(cidx,1).*pass_texit(cidx,1),numshuf,true);
    rndtenter = randsample(pass_tenter(cidx),numshuf,true);
    rndtexit = randsample(pass_texit(cidx),numshuf,true);
    crates(cvid,:) = randsample(zrate_exit(cidx,1)-zrate_enter(cidx,1),numshuf,true);
    crates_shx(cvid,:) = rndexit./rndtexit;
    crates_she(cvid,:) = rndenter./rndtenter;
end
zcrates_she = zeros(length(cvelrats),numshuf);
zcrates_shx = zeros(length(cvelrats),numshuf);
for c = 1:max(rate_enter(:,3))    
   zcrates_she(zspk_enter(:,3)==c,:) = zscore(crates_she(zspk_enter(:,3)==c,:));
   zcrates_shx(zspk_exit(:,3)==c,:) = zscore(crates_shx(zspk_exit(:,3)==c,:));      
end
crates_sh = zcrates_shx-zcrates_she;
[N,X] = hist(crates,linspace(-5,5,30));
[Nsh,Xsh] = hist(crates_sh,linspace(-5,5,30));
figure;plot(mean(Xsh,2),mean(Nsh,2),'-k');
[cish] = prctile(Nsh,[2.5 97.5],2);
hold on; shadedplot(mean(Xsh,2), cish(:,2)', cish(:,1)',[0.7 0.7 0.7],[0.5 0.5 0.5]); alpha(0.8)
hold on;plot(mean(X,2),mean(N,2),'-r');
[ci] = prctile(N,[2.5 97.5],2);
hold on; shadedplot(mean(X,2), ci(:,2)', ci(:,1)',[1 0 0],[1 0.25 0.25]); alpha(0.5)
% hold on;plot([0 td(end)],[0 0],'--k')
xlabel('exit rate / (exit rate + enter rate)');ylabel('# passes')


%% SECTION TITLE
% DESCRIPTIVE TEXT
%zscore vars by cell
zspk_enter = rate_enter;
zspk_exit = rate_exit;
zrate_enter = rate_enter;
zrate_exit = rate_exit;
zpass_tenter = pass_tenter;
zpass_texit = pass_texit;
zipass_tenter = pass_tenter;
zipass_texit = pass_texit;
zSenterdiffs = Senterdiffs;
zSexitdiffs = Sexitdiffs;
zRenterdiffs = Renterdiffs;
zRexitdiffs = Rexitdiffs;
zTenterdiffs = Tenterdiffs;
zTexitdiffs = Texitdiffs;
ziTenterdiffs = Tenterdiffs;
ziTexitdiffs = Texitdiffs;
for c = 1:max(rate_enter(:,3))    
   zspk_enter(zspk_enter(:,3)==c,1) = zscore(rate_enter(rate_enter(:,3)==c,1).*pass_tenter(pass_tenter(:,3)==c,1));
   zspk_exit(zspk_exit(:,3)==c,1) = zscore(rate_exit(rate_exit(:,3)==c,1).*pass_texit(pass_texit(:,3)==c,1));
   zrate_enter(rate_enter(:,3)==c,1) = zscore(rate_enter(rate_enter(:,3)==c,1));
   zrate_exit(rate_exit(:,3)==c,1) = zscore(rate_exit(rate_exit(:,3)==c,1));
   zpass_tenter(pass_tenter(:,3)==c,1) = zscore(pass_tenter(pass_tenter(:,3)==c,1));
   zpass_texit(pass_texit(:,3)==c,1) = zscore(pass_texit(pass_texit(:,3)==c,1));
   zipass_tenter(pass_tenter(:,3)==c,1) = zscore(1./pass_tenter(pass_tenter(:,3)==c,1));
   zipass_texit(pass_texit(:,3)==c,1) = zscore(1./pass_texit(pass_texit(:,3)==c,1));
   for p = rate_enter(rate_enter(:,3)==c,2)'       
       for s = rate_enter(rate_enter(:,3)==c & rate_enter(:,2)==p,4)'
           zSenterdiffs(Celldiffs(:,1)==c & Lapdiffs(:,1)==p & Sesdiffs(:,1)==s,1) = zspk_enter(zspk_enter(:,3)==c & zspk_enter(:,2)==p & zspk_enter(:,4)==s,1);
           zSenterdiffs(Celldiffs(:,2)==c & Lapdiffs(:,2)==p & Sesdiffs(:,2)==s,2) = zspk_enter(zspk_enter(:,3)==c & zspk_enter(:,2)==p & zspk_enter(:,4)==s,1);
           zSexitdiffs(Celldiffs(:,1)==c & Lapdiffs(:,1)==p & Sesdiffs(:,1)==s,1) = zspk_exit(zspk_exit(:,3)==c & zspk_exit(:,2)==p & zspk_exit(:,4)==s,1);
           zSexitdiffs(Celldiffs(:,2)==c & Lapdiffs(:,2)==p & Sesdiffs(:,2)==s,2) = zspk_exit(zspk_exit(:,3)==c & zspk_exit(:,2)==p & zspk_exit(:,4)==s,1);           
           zRenterdiffs(Celldiffs(:,1)==c & Lapdiffs(:,1)==p & Sesdiffs(:,1)==s,1) = zrate_enter(rate_enter(:,3)==c & rate_enter(:,2)==p & rate_enter(:,4)==s,1);
           zRenterdiffs(Celldiffs(:,2)==c & Lapdiffs(:,2)==p & Sesdiffs(:,2)==s,2) = zrate_enter(rate_enter(:,3)==c & rate_enter(:,2)==p & rate_enter(:,4)==s,1);
           zRexitdiffs(Celldiffs(:,1)==c & Lapdiffs(:,1)==p & Sesdiffs(:,1)==s,1) = zrate_exit(rate_exit(:,3)==c & rate_exit(:,2)==p & rate_exit(:,4)==s,1);
           zRexitdiffs(Celldiffs(:,2)==c & Lapdiffs(:,2)==p & Sesdiffs(:,2)==s,2) = zrate_exit(rate_exit(:,3)==c & rate_exit(:,2)==p & rate_exit(:,4)==s,1);           
           zTenterdiffs(Celldiffs(:,1)==c & Lapdiffs(:,1)==p & Sesdiffs(:,2)==s,1) = zpass_tenter(pass_tenter(:,3)==c & pass_tenter(:,2)==p & pass_tenter(:,4)==s,1);
           zTenterdiffs(Celldiffs(:,2)==c & Lapdiffs(:,2)==p & Sesdiffs(:,1)==s,2) = zpass_tenter(pass_tenter(:,3)==c & pass_tenter(:,2)==p & pass_tenter(:,4)==s,1);
           zTexitdiffs(Celldiffs(:,1)==c & Lapdiffs(:,1)==p & Sesdiffs(:,1)==s,1) = zpass_texit(pass_texit(:,3)==c & pass_texit(:,2)==p & pass_texit(:,4)==s,1);
           zTexitdiffs(Celldiffs(:,2)==c & Lapdiffs(:,2)==p & Sesdiffs(:,2)==s,2) = zpass_texit(pass_texit(:,3)==c & pass_texit(:,2)==p & pass_texit(:,4)==s,1);
           ziTenterdiffs(Celldiffs(:,1)==c & Lapdiffs(:,1)==p & Sesdiffs(:,2)==s,1) = zipass_tenter(pass_tenter(:,3)==c & pass_tenter(:,2)==p & pass_tenter(:,4)==s,1);
           ziTenterdiffs(Celldiffs(:,2)==c & Lapdiffs(:,2)==p & Sesdiffs(:,1)==s,2) = zipass_tenter(pass_tenter(:,3)==c & pass_tenter(:,2)==p & pass_tenter(:,4)==s,1);
           ziTexitdiffs(Celldiffs(:,1)==c & Lapdiffs(:,1)==p & Sesdiffs(:,1)==s,1) = zipass_texit(pass_texit(:,3)==c & pass_texit(:,2)==p & pass_texit(:,4)==s,1);
           ziTexitdiffs(Celldiffs(:,2)==c & Lapdiffs(:,2)==p & Sesdiffs(:,2)==s,2) = zipass_texit(pass_texit(:,3)==c & pass_texit(:,2)==p & pass_texit(:,4)==s,1); 
       end   
   end
end


%% SECTION TITLE
% DESCRIPTIVE TEXT
%toy poisson spiking model
rmuenter = 8.1; rventer = 4;
rmuexit = 8.7; rvexit = 5.5;
msenter = zeros(length(rate_enter(:,1)),1);
msexit = zeros(length(rate_enter(:,1)),1);
msentert1 = zeros(length(Senterdiffs(:,1)),1);
msentert2 = zeros(length(Senterdiffs(:,1)),1);
msexitt1 = zeros(length(Senterdiffs(:,1)),1);
msexitt2 = zeros(length(Senterdiffs(:,1)),1);
for c = 1:max(rate_enter(:,3))    
   Benter = regress(rate_enter(rate_enter(:,3)==c,1),[ones(size(1./pass_tenter(rate_enter(:,3)==c,1))),1./pass_tenter(rate_enter(:,3)==c,1)]);
   Bexit = regress(rate_exit(rate_enter(:,3)==c,1),[ones(size(1./pass_texit(rate_enter(:,3)==c,1))),1./pass_texit(rate_enter(:,3)==c,1)]);
   for p = rate_enter(rate_enter(:,3)==c,2)'       
       for s = rate_enter(rate_enter(:,3)==c & rate_enter(:,2)==p,4)'
            mrenter = Benter(1)+Benter(2)/pass_tenter(pass_tenter(:,2)==p&pass_tenter(:,3)==c&pass_tenter(:,4)==s,1);
            mrexit = Bexit(1)+Bexit(2)/pass_texit(pass_tenter(:,2)==p&pass_tenter(:,3)==c&pass_tenter(:,4)==s,1);
%            mrenter = mean(rate_enter(rate_enter(:,3)==c,1));
%            mrexit = mean(rate_exit(rate_enter(:,3)==c,1));
           msenter(rate_enter(:,3)==c & rate_enter(:,2)==p & rate_enter(:,4)==s) = poissrnd(mrenter*pass_tenter(zpass_tenter(:,2)==p&zpass_tenter(:,3)==c&zpass_tenter(:,4)==s,1));
           msexit(rate_enter(:,3)==c & rate_enter(:,2)==p & rate_enter(:,4)==s) = poissrnd(mrexit*pass_texit(zpass_texit(:,2)==p&zpass_texit(:,3)==c&zpass_texit(:,4)==s,1));
           msentert1(Celldiffs(:,1)==c & Lapdiffs(:,1)==p & Sesdiffs(:,1)==s) = msenter(rate_enter(:,3)==c & rate_enter(:,2)==p & rate_enter(:,4)==s);
           msentert2(Celldiffs(:,2)==c & Lapdiffs(:,2)==p & Sesdiffs(:,2)==s) = msenter(rate_enter(:,3)==c & rate_enter(:,2)==p & rate_enter(:,4)==s);
           msexitt1(Celldiffs(:,1)==c & Lapdiffs(:,1)==p & Sesdiffs(:,1)==s) = msexit(rate_enter(:,3)==c & rate_enter(:,2)==p & rate_enter(:,4)==s);
           msexitt2(Celldiffs(:,2)==c & Lapdiffs(:,2)==p & Sesdiffs(:,2)==s) = msexit(rate_enter(:,3)==c & rate_enter(:,2)==p & rate_enter(:,4)==s);
       end   
   end
end
%zscore
zmsenter = msenter;
zmsexit = msexit;
zmsentert1 = msentert1;
zmsentert2 = msentert2;
zmsexitt1 = msexitt1;
zmsexitt2 = msexitt2;
for c = 1:max(rate_enter(:,3))
    zmsenter(rate_enter(:,3)==c) = zscore(msenter(rate_enter(:,3)==c)./pass_tenter(rate_enter(:,3)==c,1));
    zmsexit(rate_enter(:,3)==c) = zscore(msexit(rate_enter(:,3)==c)./pass_texit(rate_enter(:,3)==c,1));
%     zmsenter(rate_enter(:,3)==c) = zscore(msenter(rate_enter(:,3)==c));
%     zmsexit(rate_enter(:,3)==c) = zscore(msexit(rate_enter(:,3)==c));
    for p = rate_enter(rate_enter(:,3)==c,2)'       
       for s = rate_enter(rate_enter(:,3)==c & rate_enter(:,2)==p,4)'
           zmsentert1(Celldiffs(:,1)==c & Lapdiffs(:,1)==p & Sesdiffs(:,1)==s) = zmsenter(rate_enter(:,3)==c & rate_enter(:,2)==p & rate_enter(:,4)==s);
           zmsentert2(Celldiffs(:,2)==c & Lapdiffs(:,2)==p & Sesdiffs(:,2)==s) = zmsenter(rate_enter(:,3)==c & rate_enter(:,2)==p & rate_enter(:,4)==s);
           zmsexitt1(Celldiffs(:,1)==c & Lapdiffs(:,1)==p & Sesdiffs(:,1)==s) = zmsexit(rate_enter(:,3)==c & rate_enter(:,2)==p & rate_enter(:,4)==s);
           zmsexitt2(Celldiffs(:,2)==c & Lapdiffs(:,2)==p & Sesdiffs(:,2)==s) = zmsexit(rate_enter(:,3)==c & rate_enter(:,2)==p & rate_enter(:,4)==s);
       end
    end
end

%% SECTION TITLE
% DESCRIPTIVE TEXT
%place x velocity maps, place x acceleration maps by cell
for c = 1:max(rate_enter(:,3))
    spikes = pass_spks(pass_spks(:,3)==c,1);
    spikev = pass_spkv(pass_spkv(:,3)==c,1);
    spikea = pass_spka(pass_spkv(:,3)==c,1);
    if isempty(spikes)
        continue
    end
    path = pass_path(pass_path(:,3)==c&pass_path(:,1)>=min(spikes)&pass_path(:,1)<=max(spikes),1);
    vel = pass_vel(pass_vel(:,3)==c&pass_path(:,1)>=min(spikes)&pass_path(:,1)<=max(spikes),1);
    acc = pass_acc(pass_vel(:,3)==c&pass_path(:,1)>=min(spikes)&pass_path(:,1)<=max(spikes),1);
    figure(1)
    plot(path, vel,'.k') 
    xlabel('position (rads)');ylabel('running speed (rads/sec)')
    hold on
    plot(spikes,spikev,'.r')   
    xlim([min(spikes),max(spikes)]);
    ylim([min(spikev),max(spikev)]);
    hold off
    figure(2)
    numbins = 50;
    invhl = 1/(0.1*(max(spikes)-min(spikes)));
    invh = 15;    
    map = zeros(numbins);
    lbins = linspace(min(spikes),max(spikes),numbins);
    vbins = linspace(min(spikev),max(spikev),numbins);
    for l = 1:numbins
        Kl = exp(-0.5*(spikes-lbins(l)).^2*invhl^2);
        KL = exp(-0.5*(path-lbins(l)).^2*invhl^2);
        for v = 1:numbins
            Kv = exp(-0.5*(spikev-vbins(v)).^2*invh^2);
            KV = exp(-0.5*(vel-vbins(v)).^2*invh^2);
            map(v,l) = sum(Kl.*Kv)/sum(KL.*KV);
        end
    end
    imagesc(lbins,vbins,map);axis xy;colormap hot
    xlabel('position (rads)');ylabel('running speed (rads/sec)')
%     figure(3)
%     plot(path, acc,'.k') 
%     hold on
%     plot(spikes,spikea,'.r')   
%     xlim([min(spikes),max(spikes)]);
%     ylim([min(spikea),max(spikea)]);
%     hold off
%     figure(4)
%     numbins = 50;
%     invhl = 1/(0.1*(max(spikes)-min(spikes)));
%     invh = 10;
%     map = zeros(numbins);
%     lbins = linspace(min(spikes),max(spikes),numbins);
%     abins = linspace(min(spikea),max(spikea),numbins);
%     for l = 1:numbins
%         Kl = exp(-0.5*(spikes-lbins(l)).^2*invhl^2);
%         KL = exp(-0.5*(path-lbins(l)).^2*invhl^2);
%         for a = 1:numbins
%             Ka = exp(-0.5*(spikea-abins(a)).^2*invh^2);
%             KA = exp(-0.5*(acc-abins(a)).^2*invh^2);
%             map(a,l) = sum(Kl.*Ka)/sum(KL.*KA);
%         end
%     end
%     imagesc(lbins,abins,map);axis xy;colormap hot
    pause
end


%% WPLI against coding metrics, circular correlation, etc.
% DESCRIPTIVE TEXT
bins = linspace(0,1,50);
invh = 4;
freqvec = 2:2:120;
rat = crrats(:,1);
nWPLI = numWPLI;
dWPLI = denWPLI;
PCS = zWPLIexit;
% WPLI1 = zWPLIexit-zWPLIenter;
% WPLI2 = zWPLIenter;
% WPLI = WPLI(:,pass_ct(:,4)==1);
nsh = 1000;
delta = repmat(rat,1,length(bins))-repmat(bins,length(rat),1);
% num = nWPLI*exp(-0.5*delta.*delta*invh^2);
% den = dWPLI*exp(-0.5*delta.*delta*invh^2);
% map1 = (WPLI1*exp(-0.5*delta.*delta*invh^2))./(ones(size(WPLI1))*exp(-0.5*delta.*delta*invh^2));
% map2 = (WPLI2*exp(-0.5*delta.*delta*invh^2))./(ones(size(WPLI2))*exp(-0.5*delta.*delta*invh^2));
% map = map1-map2;
map = (PCS*exp(-0.5*delta.*delta*invh^2))./(ones(size(PCS))*exp(-0.5*delta.*delta*invh^2));
% map = (num./den).^2;
mapsh = zeros(size(map,1),size(map,2),nsh);
for s = 1:nsh
    idx1 = randsample(size(PCS,2),size(PCS,2));
    %    num = nWPLI(:,idx1)*exp(-0.5*delta.*delta*invh^2);
    %    den = dWPLI(:,idx1)*exp(-0.5*delta.*delta*invh^2);
    %    sh1 = (WPLI1(:,idx1)*exp(-0.5*delta.*delta*invh^2))./(ones(size(WPLI1))*exp(-0.5*delta.*delta*invh^2));
    %    sh2 = (WPLI2(:,idx1)*exp(-0.5*delta.*delta*invh^2))./(ones(size(WPLI2))*exp(-0.5*delta.*delta*invh^2));
    %    mapsh(:,:,s) = sh1-sh2;
    mapsh(:,:,s) = (PCS(:,idx1)*exp(-0.5*delta.*delta*invh^2))./(ones(size(PCS))*exp(-0.5*delta.*delta*invh^2));
    %    mapsh(:,:,s) = (num./den).^2;
end
[ci] = prctile(mapsh,[2.5 50 97.5],3);
dmap = (map-ci(:,:,2))./std(mapsh,0,3);
% dmap(map<ci(:,:,3) & map>ci(:,:,1)) = nan;
c = max(max(dmap(freqVec>20 & freqVec<100,:)));
figure
imagesc(bins,freqvec,dmap,[-c c]);axis xy;colormap jet;colorbar
hold on; contour(bins,freqVec,dmap,[-2 2],'-k','LineWidth',2)

%% Group WPLI spectrogram
% DESCRIPTIVE TEXT
numb = 10000;
freqvec = 2:2:100;
remove = zeros(1,size(numWPLI,2));
% remove(pass_ct(:,6)==2&(pass_ct(:,5)==4|pass_ct(:,5)==5)) = 1;
% remove(pass_ct(:,6)==3&pass_ct(:,5)==2) = 1;
% remove(pass_ct(:,6)==4&(pass_ct(:,5)==1|pass_ct(:,5)==2)) = 1;
% remove(~ismember(crrats(:,3),ppc)) = 1;
nWPLI = numWPLI(:,remove==0);
dWPLI = denWPLI(:,remove==0);
idx = reshape(randsample(size(nWPLI,2),size(nWPLI,2)*numb,true),numb,size(nWPLI,2));
PCS = zeros(length(freqvec),numb);
for i = 1:numb
    PCS(:,i) = (nansum(nWPLI(:,idx(i,:)),2)./nansum(dWPLI(:,idx(i,:)),2)).^2;
end
figure;plot(freqvec,mean(PCS,2),'-k','LineWidth',1.5);
[ci] = prctile(PCS,[2.5 97.5],2);
hold on; shadedplot(freqvec, ci(:,2)', ci(:,1)',[0.7 0.7 0.7],[0.5 0.5 0.5]); alpha(0.8)
xlabel('frequency (Hz)');ylabel('WPLI')

%% Coding modes against time by session
figure; hold on;
sp = 0.5;
want1 = pass_ct(:,4)==1&pass_ct(:,7)>=0.2&pass_numspks(:,1)>=10;
want2 = pass_ct(:,4)==2&pass_ct(:,7)>=0.2&pass_numspks(:,1)>=10;
want3 = pass_ct(:,4)==3&pass_ct(:,7)>=0.2&pass_numspks(:,1)>=10;
% plot(pass_ct(want1,1)/60,crrats(want1,1),'.r');hold on
tm1 = 600/60;
plot([tm1+sp tm1+sp],[0 1],'-k');
%plot([0 tm1],[0.5 0.5],'--k')
% plot((pass_ct(want2,1))/60+tm1+2*sp,crrats(want2,1),'.r');hold on
tm2 = 600/60+tm1+2*sp;
plot([tm2+sp tm2+sp],[0 1],'-k');
%plot([tm1+2*sp tm2],[0.5 0.5],'--k')
% plot(pass_ct(want3,1)/60+tm2+2*sp,crrats(want3,1),'.r');
xlabel('time (min)');ylabel('RI');
%plot([tm2+2*sp max(pass_ct(want3,1)/60+tm2+2*sp)],[0.5 0.5],'--k')
xlim([0 600/60+tm2+2*sp]); box on

invh = 1/45;
times1 = pass_ct(want1,1);
times2 = pass_ct(want2,1);
times3 = pass_ct(want3,1);
rates1 = crrats(want1,1);
rates2 = crrats(want2,1);
rates3 = crrats(want3,1);
tbins1 = 0:0.5:600;
tbins2 = 0:0.5:600;
tbins3 = 0:0.5:600;
nsh = 1000;
RI1 = zeros(nsh,length(tbins1));
RI2 = zeros(nsh,length(tbins2));
RI3 = zeros(nsh,length(tbins3));
RI1sh = zeros(nsh,length(tbins1));
RI2sh = zeros(nsh,length(tbins2));
RI3sh = zeros(nsh,length(tbins3));
DI1 = zeros(nsh,length(tbins1));
DI2 = zeros(nsh,length(tbins2));
DI3 = zeros(nsh,length(tbins3));
DI1sh = zeros(nsh,length(tbins1));
DI2sh = zeros(nsh,length(tbins2));
DI3sh = zeros(nsh,length(tbins3));

for s = 1:nsh
   idx1 = randsample(length(times1),length(times1),true);
   idx2 = randsample(length(times2),length(times2),true);
   idx3 = randsample(length(times3),length(times3),true);   
   idx1sh = randsample(length(times1),length(times1),true);
   idx2sh = randsample(length(times2),length(times2),true);
   idx3sh = randsample(length(times3),length(times3),true);

   delta1 = repmat(times1(idx1),1,length(tbins1))-repmat(tbins1,length(times1),1);
   delta2 = repmat(times2(idx2),1,length(tbins2))-repmat(tbins2,length(times2),1);
   delta3 = repmat(times3(idx3),1,length(tbins3))-repmat(tbins3,length(times3),1);
  
   RI1(s,:) = rates1(idx1)'*exp(-0.5*delta1.*delta1*invh^2)./(ones(size(times1(idx1)))'*exp(-0.5*delta1.*delta1*invh^2));
   RI2(s,:) = rates2(idx2)'*exp(-0.5*delta2.*delta2*invh^2)./(ones(size(times2(idx2)))'*exp(-0.5*delta2.*delta2*invh^2));
   RI3(s,:) = rates3(idx3)'*exp(-0.5*delta3.*delta3*invh^2)./(ones(size(times3(idx3)))'*exp(-0.5*delta3.*delta3*invh^2));   
   RI1sh(s,:) = rates1(idx1sh)'*exp(-0.5*delta1.*delta1*invh^2)./(ones(size(times1(idx1)))'*exp(-0.5*delta1.*delta1*invh^2));
   RI2sh(s,:) = rates2(idx2sh)'*exp(-0.5*delta2.*delta2*invh^2)./(ones(size(times2(idx2)))'*exp(-0.5*delta2.*delta2*invh^2));
   RI3sh(s,:) = rates3(idx3sh)'*exp(-0.5*delta3.*delta3*invh^2)./(ones(size(times3(idx3)))'*exp(-0.5*delta3.*delta3*invh^2));  

%     fit1 = fit(times1(idx1),rates1(idx1),'smoothingspline','SmoothingParam',0.00000001);
%     fit2 = fit(times2(idx2),rates2(idx2),'smoothingspline','SmoothingParam',0.00000001);
%     fit3 = fit(times3(idx3),rates3(idx3),'smoothingspline','SmoothingParam',0.00000001);
%     fit1sh = fit(times1(idx1),rates1(idx1sh),'smoothingspline','SmoothingParam',0.00000001);
%     fit2sh = fit(times2(idx2),rates2(idx2sh),'smoothingspline','SmoothingParam',0.00000001);
%     fit3sh = fit(times3(idx3),rates3(idx3sh),'smoothingspline','SmoothingParam',0.00000001);
%     
%     RI1(s,:) = feval(fit1,tbins1)';
%     RI2(s,:) = feval(fit2,tbins2)';
%     RI3(s,:) = feval(fit3,tbins3)';
%     RI1sh(s,:) = feval(fit1sh,tbins1)';
%     RI2sh(s,:) = feval(fit2sh,tbins2)';
%     RI3sh(s,:) = feval(fit3sh,tbins3)';
%     
%     DI1(s,:) = differentiate(fit1,tbins1)';
%     DI2(s,:) = differentiate(fit2,tbins2)';
%     DI3(s,:) = differentiate(fit3,tbins3)';
%     DI1sh(s,:) = differentiate(fit1sh,tbins1)';
%     DI2sh(s,:) = differentiate(fit2sh,tbins2)';
%     DI3sh(s,:) = differentiate(fit3sh,tbins3)';
end

plot(tbins1/60,median(RI1,1),'r','LineWidth',1.5);hold on
plot(tbins2/60+tm1+2*sp,median(RI2,1),'r','LineWidth',1.5);
plot(tbins3/60+tm2+2*sp,median(RI3,1),'r','LineWidth',1.5);
plot(tbins1/60,median(RI1sh,1),'k','LineWidth',1.5);
plot(tbins2/60+tm1+2*sp,median(RI2sh,1),'k','LineWidth',1.5);
plot(tbins3/60+tm2+2*sp,median(RI3sh,1),'k','LineWidth',1.5);
[ci1] = prctile(RI1,[2.5 97.5],1);
[ci2] = prctile(RI2,[2.5 97.5],1);
[ci3] = prctile(RI3,[2.5 97.5],1);
[ci1sh] = prctile(RI1sh,[2.5 97.5],1);
[ci2sh] = prctile(RI2sh,[2.5 97.5],1);
[ci3sh] = prctile(RI3sh,[2.5 97.5],1);

shadedplot(tbins1/60, ci1sh(2,:), ci1sh(1,:),[0.7 0.7 0.7],[0.5 0.5 0.5]);alpha(0.8);hold on;
shadedplot(tbins2/60+tm1+2*sp, ci2sh(2,:), ci2sh(1,:),[0.7 0.7 0.7],[0.5 0.5 0.5]);alpha(0.8);hold on;
shadedplot(tbins3/60+tm2+2*sp, ci3sh(2,:), ci3sh(1,:),[0.7 0.7 0.7],[0.5 0.5 0.5]);alpha(0.8);hold on;
shadedplot(tbins1/60, ci1(2,:), ci1(1,:),[1 0 0],[1 0.25 0.25]); alpha(0.5);hold on;
shadedplot(tbins2/60+tm1+2*sp, ci2(2,:), ci2(1,:),[1 0 0],[1 0.25 0.25]); alpha(0.5);hold on;
shadedplot(tbins3/60+tm2+2*sp, ci3(2,:), ci3(1,:),[1 0 0],[1 0.25 0.25]); alpha(0.5);hold on;

a = gca;
% set box property to off and remove background color
set(a,'box','off','color','none')
% create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
% set original axes as active
axes(a)
% link axes in case of zooming
linkaxes([a b])

%session difference plots
figure; hold on;
sp = 0.5;
want1 = pass_ct(:,4)==1&pass_ct(:,7)>=0.2&pass_numspks(:,1)>=10;
want2 = pass_ct(:,4)==2&pass_ct(:,7)>=0.2&pass_numspks(:,1)>=10;
want3 = pass_ct(:,4)==3&pass_ct(:,7)>=0.2&pass_numspks(:,1)>=10;
% plot(pass_ct(want1,1)/60,crrats(want1,1),'.r');hold on
tm1 = 600/60;
plot([tm1+sp tm1+sp],[-0.2 0.2],'-k');
%plot([0 tm1],[0.5 0.5],'--k')
% plot((pass_ct(want2,1))/60+tm1+2*sp,crrats(want2,1),'.r');hold on
tm2 = 600/60+tm1+2*sp;
plot([tm2+sp tm2+sp],[-0.2 0.2],'-k');
%plot([tm1+2*sp tm2],[0.5 0.5],'--k')
% plot(pass_ct(want3,1)/60+tm2+2*sp,crrats(want3,1),'.r');
xlabel('time (min)');ylabel('RI');
%plot([tm2+2*sp max(pass_ct(want3,1)/60+tm2+2*sp)],[0.5 0.5],'--k')
xlim([0 600/60+tm2+2*sp]); 
ylim([-0.2 0.2]);

S1S2 = RI1 - RI2;
S1S3 = RI1 - RI3;
S2S3 = RI2 - RI3;
S1S2sh = RI1sh - RI2sh;
S1S3sh = RI1sh - RI3sh;
S2S3sh = RI2sh - RI3sh;

plot(tbins1/60,median(S1S2,1),'r','LineWidth',1.5);hold on
plot(tbins2/60+tm1+2*sp,median(S1S3,1),'r','LineWidth',1.5);
plot(tbins3/60+tm2+2*sp,median(S2S3,1),'r','LineWidth',1.5);
plot(tbins1/60,median(S1S2sh,1),'k','LineWidth',1.5);
plot(tbins2/60+tm1+2*sp,median(S1S3sh,1),'k','LineWidth',1.5);
plot(tbins3/60+tm2+2*sp,median(S2S3sh,1),'k','LineWidth',1.5);
[ci1] = prctile(S1S2,[2.5 97.5],1);
[ci2] = prctile(S1S3,[2.5 97.5],1);
[ci3] = prctile(S2S3,[2.5 97.5],1);
[ci1sh] = prctile(S1S2sh,[2.5 97.5],1);
[ci2sh] = prctile(S1S3sh,[2.5 97.5],1);
[ci3sh] = prctile(S2S3sh,[2.5 97.5],1);

shadedplot(tbins1/60, ci1sh(2,:), ci1sh(1,:),[0.7 0.7 0.7],[0.5 0.5 0.5]);alpha(0.8);hold on;
shadedplot(tbins2/60+tm1+2*sp, ci2sh(2,:), ci2sh(1,:),[0.7 0.7 0.7],[0.5 0.5 0.5]);alpha(0.8);hold on;
shadedplot(tbins3/60+tm2+2*sp, ci3sh(2,:), ci3sh(1,:),[0.7 0.7 0.7],[0.5 0.5 0.5]);alpha(0.8);hold on;
shadedplot(tbins1/60, ci1(2,:), ci1(1,:),[1 0 0],[1 0.25 0.25]); alpha(0.5);hold on;
shadedplot(tbins2/60+tm1+2*sp, ci2(2,:), ci2(1,:),[1 0 0],[1 0.25 0.25]); alpha(0.5);hold on;
shadedplot(tbins3/60+tm2+2*sp, ci3(2,:), ci3(1,:),[1 0 0],[1 0.25 0.25]); alpha(0.5);hold on;
box on

a = gca;
% set box property to off and remove background color
set(a,'box','off','color','none')
% create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
% set original axes as active
axes(a)
% link axes in case of zooming
linkaxes([a b])

% figure
% 
% plot(tbins1/60,median(DI1,1),'r','LineWidth',1.5);hold on
% plot(tbins2/60+tm1+2*sp,median(DI2,1),'r','LineWidth',1.5);
% plot(tbins3/60+tm2+2*sp,median(DI3,1),'r','LineWidth',1.5);
% plot(tbins1/60,median(DI1sh,1),'k','LineWidth',1.5);
% plot(tbins2/60+tm1+2*sp,median(DI2sh,1),'k','LineWidth',1.5);
% plot(tbins3/60+tm2+2*sp,median(DI3sh,1),'k','LineWidth',1.5);
% [ci1] = prctile(DI1,[2.5 97.5],1);
% [ci2] = prctile(DI2,[2.5 97.5],1);
% [ci3] = prctile(DI3,[2.5 97.5],1);
% [ci1sh] = prctile(DI1sh,[2.5 97.5],1);
% [ci2sh] = prctile(DI2sh,[2.5 97.5],1);
% [ci3sh] = prctile(DI3sh,[2.5 97.5],1);
% 
% shadedplot(tbins1/60, ci1sh(2,:), ci1sh(1,:),[0.7 0.7 0.7],[0.5 0.5 0.5]);alpha(0.8);hold on;
% shadedplot(tbins2/60+tm1+2*sp, ci2sh(2,:), ci2sh(1,:),[0.7 0.7 0.7],[0.5 0.5 0.5]);alpha(0.8);hold on;
% shadedplot(tbins3/60+tm2+2*sp, ci3sh(2,:), ci3sh(1,:),[0.7 0.7 0.7],[0.5 0.5 0.5]);alpha(0.8);hold on;
% shadedplot(tbins1/60, ci1(2,:), ci1(1,:),[1 0 0],[1 0.25 0.25]); alpha(0.5);hold on;
% shadedplot(tbins2/60+tm1+2*sp, ci2(2,:), ci2(1,:),[1 0 0],[1 0.25 0.25]); alpha(0.5);hold on;
% shadedplot(tbins3/60+tm2+2*sp, ci3(2,:), ci3(1,:),[1 0 0],[1 0.25 0.25]); alpha(0.5);hold on;



%% WPLI across session for backward expansion
% DESCRIPTIVE TEXT
bins = linspace(0,600,50);
invh = 1/90; sp = 60;
for b = 1:3
    rat = pass_ct(pass_ct(:,4)==b,1);
    freqvec = 2:2:120;
    freqVec = 2:2:120;
    nWPLI = numWPLI;
    dWPLI = denWPLI;
    PCS = zWPLIfull(:,pass_ct(:,4)==b);
    % WPLI1 = zWPLIexit-zWPLIenter;
    % WPLI2 = zWPLIenter;
    % WPLI = WPLI(:,pass_ct(:,4)==1);
    nsh = 1000;
    delta = repmat(rat,1,length(bins))-repmat(bins,length(rat),1);
    % num = nWPLI*exp(-0.5*delta.*delta*invh^2);
    % den = dWPLI*exp(-0.5*delta.*delta*invh^2);
    % map1 = (WPLI1*exp(-0.5*delta.*delta*invh^2))./(ones(size(WPLI1))*exp(-0.5*delta.*delta*invh^2));
    % map2 = (WPLI2*exp(-0.5*delta.*delta*invh^2))./(ones(size(WPLI2))*exp(-0.5*delta.*delta*invh^2));
    % map = map1-map2;
    map = (PCS*exp(-0.5*delta.*delta*invh^2))./(ones(size(PCS))*exp(-0.5*delta.*delta*invh^2));
    % map = (num./den).^2;
    mapsh = zeros(size(map,1),size(map,2),nsh);
    for s = 1:nsh
        idx1 = randsample(size(PCS,2),size(PCS,2));
        %    num = nWPLI(:,idx1)*exp(-0.5*delta.*delta*invh^2);
        %    den = dWPLI(:,idx1)*exp(-0.5*delta.*delta*invh^2);
        %    sh1 = (WPLI1(:,idx1)*exp(-0.5*delta.*delta*invh^2))./(ones(size(WPLI1))*exp(-0.5*delta.*delta*invh^2));
        %    sh2 = (WPLI2(:,idx1)*exp(-0.5*delta.*delta*invh^2))./(ones(size(WPLI2))*exp(-0.5*delta.*delta*invh^2));
        %    mapsh(:,:,s) = sh1-sh2;
        mapsh(:,:,s) = (PCS(:,idx1)*exp(-0.5*delta.*delta*invh^2))./(ones(size(PCS))*exp(-0.5*delta.*delta*invh^2));
        %    mapsh(:,:,s) = (num./den).^2;
    end
    [ci] = prctile(mapsh,[2.5 50 97.5],3);
    dmap = (map-ci(:,:,2))./std(mapsh,0,3);
    % dmap(map<ci(:,:,3) & map>ci(:,:,1)) = nan;
    c = max(max(dmap(freqVec>25 & freqVec<100,:)));
    imagesc(bins+(b-1)*(max(bins)+sp),freqvec,dmap,[-c c]);axis xy;colormap jet;colorbar
    hold on; contour(bins+(b-1)*(max(bins)+sp),freqVec,dmap,[-3 3],'-k','LineWidth',2)
    tm = max(bins+(b-1)*(max(bins)+sp));
    plot([tm+0.5*sp tm+0.5*sp],[0 101],'-k')
end
xlim([0 3*max(bins)+2*sp])
%% PCS across session for backward expansion
% DESCRIPTIVE TEXT
bins = linspace(0,600,50);
invh = 1/90; sp = 60;
for b = 1:3
    rat = pass_ct(pass_ct(:,4)==b&pass_ct(:,6)>3,1);
    freqvec = 2:2:100;
    freqVec = 2:2:100;    
    PCS = PCSpass(:,pass_ct(:,4)==b&pass_ct(:,6)>3);
    % WPLI1 = zWPLIexit-zWPLIenter;
    % WPLI2 = zWPLIenter;
    % WPLI = WPLI(:,pass_ct(:,4)==1);
    nsh = 1000;
    delta = repmat(rat,1,length(bins))-repmat(bins,length(rat),1);
    % num = nWPLI*exp(-0.5*delta.*delta*invh^2);
    % den = dWPLI*exp(-0.5*delta.*delta*invh^2);
    % map1 = (WPLI1*exp(-0.5*delta.*delta*invh^2))./(ones(size(WPLI1))*exp(-0.5*delta.*delta*invh^2));
    % map2 = (WPLI2*exp(-0.5*delta.*delta*invh^2))./(ones(size(WPLI2))*exp(-0.5*delta.*delta*invh^2));
    % map = map1-map2;
    map = (PCS*exp(-0.5*delta.*delta*invh^2))./(ones(size(PCS))*exp(-0.5*delta.*delta*invh^2));
    % map = (num./den).^2;
    mapsh = zeros(size(map,1),size(map,2),nsh);
    for s = 1:nsh
        idx1 = randsample(size(PCS,2),size(PCS,2));
        %    num = nWPLI(:,idx1)*exp(-0.5*delta.*delta*invh^2);
        %    den = dWPLI(:,idx1)*exp(-0.5*delta.*delta*invh^2);
        %    sh1 = (WPLI1(:,idx1)*exp(-0.5*delta.*delta*invh^2))./(ones(size(WPLI1))*exp(-0.5*delta.*delta*invh^2));
        %    sh2 = (WPLI2(:,idx1)*exp(-0.5*delta.*delta*invh^2))./(ones(size(WPLI2))*exp(-0.5*delta.*delta*invh^2));
        %    mapsh(:,:,s) = sh1-sh2;
        mapsh(:,:,s) = (PCS(:,idx1)*exp(-0.5*delta.*delta*invh^2))./(ones(size(PCS))*exp(-0.5*delta.*delta*invh^2));
        %    mapsh(:,:,s) = (num./den).^2;
    end
    [ci] = prctile(mapsh,[2.5 50 97.5],3);
    dmap = (map-ci(:,:,2))./std(mapsh,0,3);
    % dmap(map<ci(:,:,3) & map>ci(:,:,1)) = nan;
    c = max(max(dmap(freqVec>25 & freqVec<100,:)));
    imagesc(bins+(b-1)*(max(bins)+sp),freqvec,dmap,[-c c]);axis xy;colormap default;colorbar
    hold on; contour(bins+(b-1)*(max(bins)+sp),freqVec,dmap,[-3 3],'-k','LineWidth',2)
    tm = max(bins+(b-1)*(max(bins)+sp));
    plot([tm+0.5*sp tm+0.5*sp],[0 101],'-k')
    
%     if b>1
%         imagesc(bins+(b-1)*(max(bins)+sp),freqvec,map-pmap,[-0.3 0.3]);axis xy;colormap default;colorbar
%         hold on;
%         tm = max(bins+(b-1)*(max(bins)+sp));
%         plot([tm+0.5*sp tm+0.5*sp],[0 101],'-k')
%     end
%     if b==1
%         pmap = map;
%     end
end
xlim([0 3*max(bins)+2*sp])
