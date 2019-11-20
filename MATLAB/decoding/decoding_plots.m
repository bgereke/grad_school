wHPxr3 = wHPxr;
xerr3 = xerr;
wpfr3 = wpfr;
PFR3 = PFR;
wfgzmu3 = wfgzmu;
wsgzmu3 = wsgzmu;
vel3 = vel;
xhat3 = xhat;
xhatg3 = xhatg;
x3 = x;
wtfrz3 = wtfrz;
acc3 = acc;
numsp3 = numsp;
xerrg3 = xerrg;

wHPxr = [wHPxr1;wHPxr2];
xerr = [xerr1;xerr2];
wpfr = [wpfr1;wpfr2];
PFR = [PFR1;PFR2];
wfgzmu = [wfgzmu1;wfgzmu2];
wsgzmu = [wsgzmu1;wsgzmu2];
vel = [vel1;vel2];
wtfrz = [wtfrz1;wtfrz2];
x = [x1;x2];
acc = [acc1;acc2];
numsp = [numsp1;numsp2];
xerrg = [xerrg1;xerrg2];

wHPxr = [wHPxr1;wHPxr2;wHPxr3];
xerr = [xerr1;xerr2;xerr3];
xerrg = [xerrg1;xerrg2;xerrg3];
wpfr = [wpfr1;wpfr2;wpfr3];
PFR = [PFR1;PFR2;PFR3];
wfgzmu = [wfgzmu1;wfgzmu2;wfgzmu3];
wsgzmu = [wsgzmu1;wsgzmu2;wsgzmu3];
vel = [vel1;vel2;vel3];
wtfrz = [wtfrz1;wtfrz2;wtfrz3];
x = [x1;x2;x3];
acc = [acc1;acc2;acc3];
numsp = [numsp1;numsp2;numsp3];
xerrg = [xerrg1;xerrg2;xerrg3];
xhat = [xhat1;xhat2;xhat3];
xhatg = [xhatg1;xhatg2;xhatg3];


H = []; for i=1:size(wHPxr,1), H = [H;wHPxr{i}]; end
E = []; for i=1:size(wHPxr,1), E = [E;xerr{i}]; end
Eg = []; for i=1:size(wHPxr,1), Eg = [Eg;xerrg{i}]; end
V = []; for i=1:size(wHPxr,1), V = [V;vel(i)]; end
X = []; for i=1:size(wHPxr,1), X = [X;xhat{i}]; end
FG = []; for i=1:size(wHPxr,1), FG = [FG;wfgzmu{i}]; end
SG = []; for i=1:size(wHPxr,1), SG = [SG;wsgzmu{i}]; end
FGz = (FG - nanmean(FG))/(nanvar(FG));
SGz = (SG - nanmean(SG))/(nanvar(SG));
Hz = (H - nanmean(H))/(nanvar(H));
Ez = (abs(E(:,2)) - nanmean(abs(E(:,2))))/(nanvar(abs(E(:,2))));
Vz = (abs(V) - nanmean(abs(V)))/(nanvar(abs(V)));


Ehp = E(abs(V)>10,:);
Hhp = H(abs(V)>10);
Vhp = V(abs(V)>10);
FGhp = FG(abs(V)>10);
SGhp = SG(abs(V)>10);
wpfrhp = wpfr{abs(V)>10};

freqvec = freqVec;

%Low-High Entropy
P1 = 0; n = 0; m = 0;
Hn = -H;
for i=1:size(H,1)
    if H(i) < 2
        n = n+1; m = m + Hn(i);
        P1(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i};
    end
end
P1 = nansum(P1,3)/n;n1=n
P2 = 0; n = 0; m = 0;
Hn = H;
for i=1:size(H,1)
    if H(i) > 2
        n = n+1; m = m + Hn(i);
        P2(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i};
    end
end
P2 = nansum(P2,3)/n;n2=n
P = P1-P2;
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P1,[-0.4 0.4]); axis xy; colorbar
title('Low Entropy');
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P2,[-0.4 0.4]); axis xy; colorbar
title('High Entropy');
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P,[-0.2 0.2]); axis xy; colorbar
title('Low-High Entropy');
xlabel('theta phase'); ylabel('Frequency (Hz)');

%Retro-Prospective
% P1 = 0; n = 0; for i=1:size(wHPxr,1), if xerr{i}(2)<0,n = n+1; P1(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i};  end, end
% P1 = nansum(P1,3)/n;n1=n
% P2 = 0; n = 0; for i=1:size(wHPxr,1), if xerr{i}(2)>0,n = n+1; P2(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i};  end, end
% P2 = nansum(P2,3)/n;n2=n
% P = P1-P2;
% figure
% imagesc(-pi:10/360*2*pi:pi,freqvec,P1,[-0.4 0.4]); axis xy; colorbar
% title('Retrospective');
% figure
% imagesc(-pi:10/360*2*pi:pi,freqvec,P2,[-0.4 0.4]); axis xy; colorbar
% title('Prospective');
% figure
% imagesc(-pi:10/360*2*pi:pi,freqvec,P,[-0.2 0.2]); axis xy; colorbar
% title('Retro-Prospective');
% xlabel('theta phase'); ylabel('Frequency (Hz)');

%Retro-Prospective (proportional to error)
P1 = 0; n = 0; m = 0;
for i=1:size(H,1) 
    if E(i,2)<-25 
        n = n+1; m = m + log2(abs(E(i,2)));
        P1(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i};  
    end
end
P1 = nansum(P1,3)/n;n1=n

P2 = 0; n = 0; m =0; 
for i=1:size(H,1)
    if E(i,2)>25
        n = n+1; m = m+ log(abs(E(i,2)));
        P2(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i};  
    end
end
P2 = nansum(P2,3)/n;n2=n
P = P1-P2;
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P1,[-0.4 0.4]); axis xy; colorbar
title('x* - x < -75');xlabel('theta phase'); ylabel('Frequency (Hz)');
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P2,[-0.4 0.4]); axis xy; colorbar
title('x* - x > 75');xlabel('theta phase'); ylabel('Frequency (Hz)');
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P,[-0.4 0.4]); axis xy; colorbar
title('Retrospective - Prospective');
xlabel('theta phase'); ylabel('Frequency (Hz)');

%Retro-Prospective + Low Entropy
P1 = 0; n = 0;
for i=1:size(H,1) 
    if E(i,2)<0 & H(i) < nanmedian(H),
        n = n+1;
        P1(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i};  
    end
end
P1 = nansum(P1,3)/n;n1=n

P2 = 0; n = 0;
for i=1:size(H,1)
    if E(i,2)>0 & H(i) < nanmedian(H)
        n = n+1;
        P2(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i};  
    end
end
P2 = nansum(P2,3)/n;n2=n
P = P1-P2;
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P1,[-0.4 0.4]); axis xy; colorbar
title('Retrospective');
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P2,[-0.4 0.4]); axis xy; colorbar
title('Prospective');
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P,[-0.2 0.2]); axis xy; colorbar
title('Retro-Prospective + Low Entropy');
xlabel('theta phase'); ylabel('Frequency (Hz)');

%Retro-Prospective + Low Entropy (proportional to error)
P1 = 0; n = 0; m = 0;
for i=1:size(Hhp,1) 
    if xerr{i}(2)<0 & H(i) < nanmedian(H),
        n = n+1; m = m + abs(xerr{i}(2)^2);
        P1(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i}*abs(xerr{i}(2)^2);  
    end
end
P1 = nansum(P1,3)/m;n1=n

P2 = 0; n = 0; m =0; 
for i=1:size(Hhp,1)
    if xerr{i}(2)>0 & H(i) < nanmedian(H)
        n = n+1; m = m+ abs(xerr{i}(2)^2);
        P2(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i}*abs(xerr{i}(2)^2);  
    end
end
P2 = nansum(P2,3)/m;n2=n
P = P1-P2;
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P1,[-0.4 0.4]); axis xy; colorbar
title('Retrospective');
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P2,[-0.4 0.4]); axis xy; colorbar
title('Prospective');
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P,[-0.2 0.2]); axis xy; colorbar
title('Retro-Prospective + Low Entropy (proportion)');
xlabel('theta phase'); ylabel('Frequency (Hz)');

%Retro-Prospective + Low Entropy (inversely proportional to error)
P1 = 0; n = 0; m = 0;
for i=1:size(Hhp,1) 
    if xerr{i}(2)<0 & H(i) < nanmedian(H),
        n = n+1; m = m + 1/abs(xerr{i}(2)^2);
        P1(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i}/abs(xerr{i}(2)^2);  
    end
end
P1 = nansum(P1,3)/m;n1=n

P2 = 0; n = 0; m =0; 
for i=1:size(Hhp,1)
    if xerr{i}(2)>0 & H(i) < nanmedian(H)
        n = n+1; m = m+ 1/abs(xerr{i}(2)^2);
        P2(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i}/abs(xerr{i}(2)^2);  
    end
end
P2 = nansum(P2,3)/m;n2=n
P = P1-P2;
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P1,[-0.4 0.4]); axis xy; colorbar
title('Retrospective');
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P2,[-0.4 0.4]); axis xy; colorbar
title('Prospective');
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P,[-0.2 0.2]); axis xy; colorbar
title('Retro-Prospective + Low Entropy (inverse proportion)');
xlabel('theta phase'); ylabel('Frequency (Hz)');

%Retro-Prospective + High Entropy
P1 = 0; n = 0; for i=1:size(Hhp,1), if xerr{i}(2)<0 & H(i) > nanmedian(H),n = n+1; P1(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i};  end, end
P1 = nansum(P1,3)/n;n1=n
P2 = 0; n = 0; for i=1:size(Hhp,1), if xerr{i}(2)>0 & H(i) > nanmedian(H),n = n+1; P2(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i};  end, end
P2 = nansum(P2,3)/n;n2=n
P = P1-P2;
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P1,[-0.4 0.4]); axis xy; colorbar
title('Retrospective');
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P2,[-0.4 0.4]); axis xy; colorbar
title('Prospective');
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P,[-0.2 0.2]); axis xy; colorbar
title('Retro-Prospective + High Entropy');
xlabel('theta phase'); ylabel('Frequency (Hz)');

%Small-Large Error 
P1 = 0; n = 0; m = 0;
for i=1:size(Hhp,1)
    if ~isnan(xerr{i}(2))
        n = n+1; m = m + 1/abs(xerr{i}(2)^2);
        P1(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i}/abs(xerr{i}(2)^2);
    end
end
P1 = nansum(P1,3)/m;n1=n
P2 = 0; n = 0; m = 0;
for i=1:size(Hhp,1)
    if ~isnan(xerr{i}(2))
        n = n+1; m = m + abs(xerr{i}(2)^2);
        P2(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i}*abs(xerr{i}(2)^2);  
    end
end
P2 = nansum(P2,3)/m;n2=n
P = P1-P2;figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P1,[-0.4 0.4]); axis xy; colorbar
title('Small');
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P2,[-0.4 0.4]); axis xy; colorbar
title('Large');
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P,[-0.2 0.2]); axis xy; colorbar
title('Small-Large Error');
xlabel('theta phase'); ylabel('Frequency (Hz)');

%High Slow Gamma
P1 = 0; n = 0; for i=1:size(Hhp,1), if SG(i)>2,n = n+1; P1(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i};  end, end
P1 = nansum(P1,3)/n;n1=n
P2 = 0; n = 0; for i=1:size(Hhp,1), if SG(i)<1,n = n+1; P2(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i};  end, end
P2 = nansum(P2,3)/n;n2=n
P = P1-P2;figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P1,[-0.5 0.5]); axis xy; colorbar
title('Large Slow Gamma');
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P2,[-0.5 0.5]); axis xy; colorbar
title('Small Slow Gamma');
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P,[-0.5 0.5]); axis xy; colorbar
title('Large-Small Slow Gamma');
xlabel('theta phase'); ylabel('Frequency (Hz)');

%High Fast Gamma
P1 = 0; n = 0; for i=1:size(Hhp,1), if FG(i)>2,n = n+1; P1(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i};  end, end
P1 = nansum(P1,3)/n;n1=n
P2 = 0; n = 0; for i=1:size(Hhp,1), if FG(i)<2,n = n+1; P2(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i};  end, end
P2 = nansum(P2,3)/n;n2=n
P = P1-P2;figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P1,[-0.5 0.5]); axis xy; colorbar
title('Large Fast Gamma');
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P2,[-0.5 0.5]); axis xy; colorbar
title('Small Fast Gamma');
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P,[-0.5 0.5]); axis xy; colorbar
title('Large-Small Fast Gamma');
xlabel('theta phase'); ylabel('Frequency (Hz)');

%Fast-Slow Gamma
P1 = 0; n = 0; for i=1:size(Hhp,1), if FG(i)>2,n = n+1; P1(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i};  end, end
P1 = nansum(P1,3)/n;n1=n
P2 = 0; n = 0; for i=1:size(Hhp,1), if SG(i)>2,n = n+1; P2(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i};  end, end
P2 = nansum(P2,3)/n;n2=n
P = P1-P2;figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P1,[-0.5 0.5]); axis xy; colorbar
title('Large Fast Gamma');
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P2,[-0.5 0.5]); axis xy; colorbar
title('Large Slow Gamma');
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P,[-0.5 0.5]); axis xy; colorbar
title('Fast-Slow Gamma');
xlabel('theta phase'); ylabel('Frequency (Hz)');

%High-Low Velocity
P1 = 0; n = 0; m = 0;
Vn = log2(abs(V));
for i=1:size(H,1)
    if Vn(i) > 3
        n = n+1; m = m + Vn(i);
        P1(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i};
    end
end
P1 = nansum(P1,3)/n;n1=n
P2 = 0; n = 0; m = 0;
%Vn = -abs(Vhp) + max(abs(Vhp)) + min(abs(Vhp));
for i=1:size(H,1)
    if Vn(i) < 3 
        n = n+1; m = m + Vn(i);
        P2(1:size(wpfr{i},1),1:size(wpfr{i},2),n) = wpfr{i};
    end
end
P2 = nansum(P2,3)/n;n2=n
P = P1-P2;
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P1,[-0.4 0.4]); axis xy; colorbar
title('Log-Velocity > 3');xlabel('theta phase'); ylabel('Frequency (Hz)');
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P2,[-0.4 0.4]); axis xy; colorbar
title('Log-Velocity < 3');xlabel('theta phase'); ylabel('Frequency (Hz)');
figure
imagesc(-pi:10/360*2*pi:pi,freqvec,P,[-0.3 0.3]); axis xy; colorbar
title('Moving-Not Moving Log-Velocity');
xlabel('theta phase'); ylabel('Frequency (Hz)');

%Negative Error Triggered Covariance
X = zeros(size(wpfr{1},1),size(wpfr{1},2),size(wpfr,1));
y = E(:,2) < 0;
for i = 1:size(wpfr,1)
    X(:,:,i) = wpfr{i}; 
end
[STA,STC,UC,SC,VC,U,S,V] = STA_STC(X,y);

%Velocity Frequency Representation
freqvec = freqVec; V = vel;
velvec = 1:1:20;
velbins = velvec;
vfr = zeros(length(freqvec),length(velbins)-1);
for vb=1:length(velbins)-1    
           vidx = find(abs(V) >= velbins(vb) & abs(V) < velbins(vb+1)); 
           wtfrz_b = []; wtd_b = [];
           for i = 1:length(vidx)      
                 wtfrz_b = [wtfrz_b wtfrz{vidx(i)}];
                 %wtd_b = [wtd_b wthetadelta{vidx(i)}];
                 %vfr(:,vb) = vfr(:,vb) + mean(wtfrz{vidx(i)},2)/length(vidx);                             
           end
           vfr(:,vb) = mean(wtfrz_b,2);
           %vfr(:,vb) = mean(wtfrz_b(:,wtd_b>2),2);
end
figure
uimagesc(velbins(1:end-1),freqvec,vfr,[-0.3,0.3]);axis xy; colorbar
xlabel('Velocity');ylabel('Frequency (Hz)');title('Velocity Frequency Representation')
set(gca,'xtick',velbins(1:3:end))
set(gca,'XTickLabel',velvec(1:3:end))
a = freqvec;
set(gca,'ytick',a(1:10:end))
set(gca,'YTickLabel',a(1:10:end))

%Acceleration Frequency Representation
accvec = -20:2:20;
accbins = accvec;
afr = zeros(length(freqVec),length(accbins)-1);
for ab=1:length(accbins)-1    
           aidx = find(acc >= accbins(ab) & acc < accbins(ab+1)); 
           wtfrz_b = [];wtd_b = [];
           for i = 1:length(aidx)     
               wtfrz_b = [wtfrz_b wtfrz{aidx(i)}];
               %wtd_b = [wtd_b wthetadelta{aidx(i)}];
                %afr(:,ab) = afr(:,ab) + mean(wtfrz{aidx(i)},2)/length(aidx);                             
           end
           afr(:,ab) = mean(wtfrz_b,2);
          % afr(:,ab) = mean(wtfrz_b(:,wtd_b>2),2);
end
figure
imagesc(accbins(1:end-1),freqVec,afr,[-0.3 0.3]);axis xy; colorbar
xlabel('Acceleration');ylabel('Frequency (Hz)');title('Acceleration Frequency Representation')
set(gca,'xtick',log2(velbins(1:3:end)))
set(gca,'XTickLabel',velvec(1:3:end))
a = freqVec;
set(gca,'ytick',a(1:10:end))
set(gca,'YTickLabel',a(1:10:end))

%Velocity Spectral Correlation
freqvec = freqVec; V = vel;
velvec = 1:1:25;
velbins = velvec;
vfr = zeros(length(freqvec),length(velbins)-1);
for vb=1:length(velbins)-1    
           vidx = find(abs(V) >= velbins(vb) & abs(V) < velbins(vb+1)); 
           wtfrz_b1 = []; wtfrz_b2 = []; wtd_b = [];
           for i = 1:length(vidx)      
                 wtfrz_b1 = [wtfrz_b1 wtfrz1{vidx(i)}];
                 wtfrz_b2 = [wtfrz_b2 wtfrz2{vidx(i)}];
                 %wtd_b = [wtd_b wthetadelta{vidx(i)}];
                 %vfr(:,vb) = vfr(:,vb) + mean(wtfrz{vidx(i)},2)/length(vidx);                             
           end           
           R = zeros(size(freqvec));
           for i = 1:length(freqvec)
               r = corrcoef(wtfrz_b1(i,:),wtfrz_b2(i,:));
               R(i) = r(1,2);
           end
           vfr(:,vb) = R;
           %vfr(:,vb) = mean(wtfrz_b(:,wtd_b>2),2);
end
figure
uimagesc(velbins(1:end-1),freqvec,vfr);axis xy; colorbar
xlabel('Velocity');ylabel('Frequency (Hz)');title('Velocity Spectral Correlation')
set(gca,'xtick',velbins(1:3:end))
set(gca,'XTickLabel',velvec(1:3:end))
a = freqvec;
set(gca,'ytick',a(1:10:end))
set(gca,'YTickLabel',a(1:10:end))

%Acceleration Spectral Correlation
accvec = -20:2:20;
accbins = accvec;
afr = zeros(length(freqVec),length(accbins)-1);
for ab=1:length(accbins)-1    
           aidx = find(acc >= accbins(ab) & acc < accbins(ab+1)); 
           wtfrz_b1 = [];wtfrz_b2 = [];wtd_b = [];
           for i = 1:length(aidx)     
                wtfrz_b1 = [wtfrz_b1 wtfrz1{aidx(i)}];
                wtfrz_b2 = [wtfrz_b2 wtfrz2{aidx(i)}];
               %wtd_b = [wtd_b wthetadelta{aidx(i)}];
                %afr(:,ab) = afr(:,ab) + mean(wtfrz{aidx(i)},2)/length(aidx);                             
           end
            R = zeros(size(freqvec));
           for i = 1:length(freqvec)
               r = corrcoef(wtfrz_b1(i,:),wtfrz_b2(i,:));
               R(i) = r(1,2);
           end
           afr(:,ab) = R;
          % afr(:,ab) = mean(wtfrz_b(:,wtd_b>2),2);
end
figure
imagesc(accbins(1:end-1),freqVec,afr);axis xy; colorbar
xlabel('Acceleration');ylabel('Frequency (Hz)');title('Acceleration Spectral Correlation')
set(gca,'xtick',log2(velbins(1:3:end)))
set(gca,'XTickLabel',velvec(1:3:end))
a = freqVec;
set(gca,'ytick',a(1:10:end))
set(gca,'YTickLabel',a(1:10:end))

%Spike Frequency Representation
spvec = 2:23;
spbins = spvec;
spfr = zeros(length(freqVec),length(spbins)-1);
for spb=1:length(spbins)-1    
           spidx = find(numsp >= spbins(spb) & numsp < spbins(spb+1));           
           for i = 1:length(spidx)               
                 spfr(:,spb) = spfr(:,spb) + mean(wtfrz{spidx(i)},2)/length(spidx);                             
           end
end
figure
imagesc(spbins(1:end-1),freqVec,spfr,[-0.3 0.3]);axis xy; colorbar
xlabel('# Spikes');ylabel('Frequency (Hz)');title('Spike Frequency Representation')
set(gca,'xtick',log2(velbins(1:3:end)))
set(gca,'XTickLabel',velvec(1:3:end))
a = freqVec;
set(gca,'ytick',a(1:10:end))
set(gca,'YTickLabel',a(1:10:end))

%Error-Velocity Scatter Plot
figure
scatter(log2(abs(V)),log2(abs(E(:,2))))
xlabel('Log-Velocity');ylabel('Log-Error');title('Absolute Error vs Velocity')

%Error Frequency Representation
evec = min(E(abs(E(:,2))<45,2)):3:max(E(abs(E(:,2))<45,2));
errbins = evec;
efr = zeros(length(freqvec),length(errbins)-1);
for eb=1:size(errbins,2)-1    
           eidx = find(E(:,2) >= errbins(eb) & E(:,2) < errbins(eb+1));
           for i = 1:length(eidx)
                    efr(:,eb) = efr(:,eb) + mean(wtfrz{eidx(i)},2)/length(eidx);               
           end
 
end
figure
uimagesc(errbins(1:end-1),freqvec,efr,[-0.2 0.2]);axis xy; colorbar
xlabel('x* - x');ylabel('frequency (Hz)');title('Error Frequency Representation')

%Absolute Error Frequency Representation
evec = min(abs(E(:,2))):3:max(E(abs(E(:,2))<45,2));
errbins = evec;
efr = zeros(length(freqvec),length(errbins)-1);
for eb=1:size(errbins,2)-1    
           eidx = find(abs(E(:,2)) >= errbins(eb) & abs(E(:,2)) < errbins(eb+1));
           for i = 1:length(eidx)
                    efr(:,eb) = efr(:,eb) + mean(wtfrz{eidx(i)},2)/length(eidx);               
           end
 
end
figure
uimagesc(errbins(1:end-1),freqvec,efr,[-0.1 0.1]);axis xy; colorbar
xlabel('|x* - x|');ylabel('frequency (Hz)');title('Absolute Error Frequency Representation')

%Space Frequency Representation
xvec = xbins;
xbins = xvec;
xfr = zeros(length(freqvec),length(xbins));
for xb=1:length(xbins)   
           xidx = find(x == xbins(xb));
           for i = 1:length(xidx)
                    xfr(:,xb) = xfr(:,xb) + mean(wtfrz{xidx(i)},2)/length(xidx);               
           end
 
end
figure
uimagesc(xbins(1:end),freqvec,xfr,[-0.25 0.25]);axis xy; colorbar
xlabel('position (cm)');ylabel('frequency (Hz)');title('Space Frequency Representation')

%Entropy Velocity Scatter
figure
scatter(abs(V),H)
xlabel('Velocity');ylabel('H');title('Entropy vs Velocity')

%Entropy Frequency Representation
hvec = 0.1:0.25:5.5;
hbins = hvec;
hfr = zeros(length(freqvec),length(hbins)-1);
for hb=1:size(hbins,2)-1    
           hidx = find(H >= hbins(hb) & H < hbins(hb+1));
           for i = 1:length(hidx)
                    hfr(:,hb) = hfr(:,hb) + mean(wtfrz{hidx(i)},2)/length(hidx);               
           end
 
end
figure
uimagesc(log2(hbins(1:end-1)),log2(freqvec),hfr,[-0.2 0.2]);axis xy; colorbar
xlabel('H');ylabel('frequency (Hz)');title('Entropy Frequency Representation')
set(gca,'xtick',log2(hbins(1:3:end)))
set(gca,'XTickLabel',hvec(1:3:end))
a = freqvec;
set(gca,'ytick',log2(a(1:10:end)))
set(gca,'YTickLabel',a(1:10:end))

%Entropy Frequency Representation
hbins = 0.1:0.15:5.5;
hfr = zeros(length(freqvec),length(hbins)-1);
for hb=1:size(hbins,2)-1    
           hidx = find(H >= hbins(hb) & H < hbins(hb+1) & abs(vel));
           for i = 1:length(hidx)
                    hfr(:,hb) = hfr(:,hb) + mean(wtfrz{hidx(i)},2)/length(hidx);               
           end
 
end
figure
uimagesc(hbins(1:end-1),freqvec,hfr,[-0.3 0.3]);axis xy; colorbar
xlabel('H');ylabel('frequency (Hz)');title('Entropy Frequency Representation')
set(gca,'xtick',hbins(1:3:end))
set(gca,'XTickLabel',hbins(1:3:end))
a = freqvec;
set(gca,'ytick',a(1:10:end))
set(gca,'YTickLabel',a(1:10:end))

%Error vs Entropy
figure
scatter(H,abs(E(:,2)))
xlabel('H');ylabel('|x* - x|');title('Error vs Entropy')

%Relating Error to Velocity
[vt vc] = hist(abs(V),100);
VC = zeros(length(V),1);
for i = 1:length(V)
    [vmin vcidx]  = min((vc - abs(V(i))).^2);
    VC(i) = vt(vcidx);
end
Ez = (abs(E(:,2))-nanmean(abs(E(:,2))))/(nanvar(abs(E(:,2))));
VCz = zscore(VC);

scatter(abs(V),abs(Ez./VCz));

