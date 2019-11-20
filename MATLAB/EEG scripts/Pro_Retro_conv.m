%Spatial Convolution of mode-triggered fields
cd('C:\Data\Pro_Retro_data\')
load('data.mat')

%create directory to store files
dd = cd; %data directory
dirInfo = dir(cd);
found1 = 0;
found2 = 0;
found3 = 0;
found4 = 0;
found5 = 0;
found6 = 0;
for ss=1:size(dirInfo,1)
    if dirInfo(ss).isdir
        if strcmp(dirInfo(ss).name,strcat('Convolution_rmaps','\'))
            found1 = 1;
        end
        if strcmp(dirInfo(ss).name,strcat('Convolution_pmaps','\'))
            found2 = 1;
        end     
        if strcmp(dirInfo(ss).name,strcat('Convolution_amaps_first','\'))
            found3 = 1;
        end
        if strcmp(dirInfo(ss).name,strcat('Convolution_rinfo','\'))
            found4 = 1;
        end    
        if strcmp(dirInfo(ss).name,strcat('Convolution_pinfo','\'))
            found5 = 1;
        end
        if strcmp(dirInfo(ss).name,strcat('Convolution_ainfo_first','\'))
            found6 = 1;
        end  
    end
end
if found1==0
    mkdir(dd,strcat('\Convolution_rmaps','\'));
end
if found2==0
    mkdir(dd,strcat('\Convolution_pmaps','\'));
end
if found3==0
    mkdir(dd,strcat('\Convolution_amaps_first','\'));
end
if found4==0
    mkdir(dd,strcat('\Convolution_rinfo','\'));
end
if found5==0
    mkdir(dd,strcat('\Convolution_pinfo','\'));
end
if found6==0
    mkdir(dd,strcat('\Convolution_ainfo_first','\'));
end

cells = unique(cmodes(:,3));
numcells = length(cells);
tshift = -3:0.02:3;      %time shifts for convolution
D = 100;    %cm
binsize = 5; %cm
% mapAxis = linspace(-pi,pi,pi*D/binsize);
mapAxis = linspace(-2,2,150);
dx = mapAxis(2)-mapAxis(1);
kappa = 100; % Smoothing factor when calculating the ratemap
pcut = 30;  %number of passes required to compute rate maps
acut = 25;
rpeak = [];
rpeakt = [];
rskew = [];
ppeak = [];
ppeakt = [];
pskew = [];
apeak = [];
apeakt = [];
askew = [];

for nc = 1:numcells
    
    convr = zeros(length(mapAxis),length(tshift));  %mode-triggered convolution maps
    convp = convr; conva = convr; 
    rinfo = zeros(length(tshift),1);     %mode triggered spatial information
    pinfo = rinfo; ainfo = rinfo;
    
    numr = sum(cmodes(:,3)==cells(nc)&cmodes(:,1)==1);    %number of mode-triggered passes
    nump = sum(cmodes(:,3)==cells(nc)&cmodes(:,1)==-1);
    numa = sum(cmodes(:,3)==cells(nc)&isnan(cmodes(:,1)));
    
    if numr >= pcut 
        
        rpasses = find(cmodes(:,3)==cells(nc)&cmodes(:,1)==1);        
        rts = [];   %mode-triggered spike times        
        rpath = [];     %mode triggered animal trajectories        
        rt = [];        %mode-triggered tracking times
        
        for rp = rpasses'
            rRows = ismember(pass_ts(:,2:6),cmodes(rp,2:6),'rows');
            rts = [rts;pass_ts(rRows,1)];
            rRows = ismember(pass_path(:,2:6),cmodes(rp,2:6),'rows');
            rpath = [rpath;pass_path(rRows,1)];
            rt = [rt;pass_t(rRows,1)];
        end
        
        for ts = 1:length(tshift)            
            rts_sh = rts + tshift(ts);
            [rspk_phase] = spikePos(rts_sh,rpath,rt);
            [convr(:,ts),rtmap,~] = circle_map(rspk_phase,rpath',rt',kappa,mapAxis); %conv = rate map, rtmap = time map 
            rinfo(ts) = (rtmap/sum(rtmap).*convr(:,ts)/(convr(:,ts)'*rtmap/sum(rtmap)))'*log2(convr(:,ts)/(convr(:,ts)'*rtmap/sum(rtmap))); %spatial information (bit/spike)
        end
        
        [p,pidx] = max(rinfo);
        rpeak = [rpeak;p];
        rpeakt = [rpeakt;tshift(pidx)];
        rskew = [rskew;myskewness(tshift,0,[],rinfo)];       
        
        figure(1)
        imagesc(tshift,mapAxis,convr);axis xy;colorbar;colormap hot
        xlabel('times shift (sec)');ylabel('track position (rads)');
        title(strcat('Convolution of Retrospective Laps: Cell ',num2str(cells(nc))));
        hold on;plot([0 0],[min(mapAxis) max(mapAxis)],'b');plot([min(tshift) max(tshift)],[0 0],'b');hold off
        figImage = sprintf('%s%s%s%s',cd,strcat('\Convolution_rmaps','\'),strcat('cell_',num2str(cells(nc))),'.fig');
        bmpImage = sprintf('%s%s%s%s',cd,strcat('\Convolution_rmaps','\'),strcat('cell_',num2str(cells(nc))),'.bmp');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');       
        plot(tshift,rinfo,'r');hold on;
        xlabel('time shift (sec)');
        ylabel('spatial information (bits/spike)');
        title(strcat('Spatial Information Against Spike Shifts: Cell ',num2str(cells(nc))));
        legend('retrospective laps');hold off;
        figImage = sprintf('%s%s%s%s',cd,strcat('\Convolution_rinfo','\'),strcat('cell_',num2str(cells(nc))),'.fig');
        bmpImage = sprintf('%s%s%s%s',cd,strcat('\Convolution_rinfo','\'),strcat('cell_',num2str(cells(nc))),'.bmp');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');
        close all
    end
    if nump >= pcut 
         ppasses = find(cmodes(:,3)==cells(nc)&cmodes(:,1)==-1);
         pts = [];      %mode triggered spike times
         ppath = [];    %mode triggered animal trajectories 
         pt = [];       %mode-triggered tracking times  
         
         for pp = ppasses'
             pRows = ismember(pass_ts(:,2:6),cmodes(pp,2:6),'rows');
             pts = [pts;pass_ts(pRows,1)];
             pRows = ismember(pass_path(:,2:6),cmodes(pp,2:6),'rows');
             ppath = [ppath;pass_path(pRows,1)];
             pt = [pt;pass_t(pRows,1)];
         end
         for ts = 1:length(tshift)
             pts_sh = pts + tshift(ts);
             [pspk_phase] = spikePos(pts_sh,ppath,pt);
             [convp(:,ts),ptmap,~] = circle_map(pspk_phase,ppath',pt',kappa,mapAxis);
             pinfo(ts) = (ptmap/sum(ptmap).*convp(:,ts)/(convp(:,ts)'*ptmap/sum(ptmap)))'*log2(convp(:,ts)/(convp(:,ts)'*ptmap/sum(ptmap))); %spatial information (bit/spike)
         end
         
        [p,pidx] = max(pinfo);
        ppeak = [ppeak;p];
        ppeakt = [ppeakt;tshift(pidx)];
        pskew = [pskew;myskewness(tshift,0,[],pinfo)];       
        
        figure(1)
        imagesc(tshift,mapAxis,convp);axis xy;colorbar;colormap hot
        xlabel('times shift (sec)');ylabel('track position (rads)');
        title(strcat('Convolution of Prospective Laps: Cell ',num2str(cells(nc))));
        hold on;plot([0 0],[min(mapAxis) max(mapAxis)],'b');plot([min(tshift) max(tshift)],[0 0],'b');hold off
        figImage = sprintf('%s%s%s%s',cd,strcat('\Convolution_pmaps','\'),strcat('cell_',num2str(cells(nc))),'.fig');
        bmpImage = sprintf('%s%s%s%s',cd,strcat('\Convolution_pmaps','\'),strcat('cell_',num2str(cells(nc))),'.bmp');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');       
        plot(tshift,pinfo,'b');hold on;
        xlabel('time shift (sec)');
        ylabel('spatial information (bits/spike)');
        title(strcat('Spatial Information Against Spike Shifts: Cell ',num2str(cells(nc))));
        legend('prospective laps');hold off;
        figImage = sprintf('%s%s%s%s',cd,strcat('\Convolution_pinfo','\'),strcat('cell_',num2str(cells(nc))),'.fig');
        bmpImage = sprintf('%s%s%s%s',cd,strcat('\Convolution_pinfo','\'),strcat('cell_',num2str(cells(nc))),'.bmp');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');
        close all
    end
    if numa+numr+nump >= acut 
        apasses = find(cmodes(:,3)==cells(nc));%&isnan(cmodes(:,1)));
        ats = [];   %mode triggered spike times
        apath = []; %mode triggered animal trajectories 
        at = [];    %mode-triggered tracking times
        
        for ap = apasses'
            aRows = ismember(pass_ts(:,2:6),cmodes(ap,2:6),'rows');
            idx = find(aRows,1,'first'); %first spike only
            ats = [ats;pass_ts(idx,1)];
            aRows = ismember(pass_path(:,2:6),cmodes(ap,2:6),'rows');
            apath = [apath;pass_path(aRows,1)];
            at = [at;pass_t(aRows,1)];
        end 
        for ts = 1:length(tshift)            
            ats_sh = ats + tshift(ts);
            [aspk_phase] = spikePos(ats_sh,apath,at); 
%             scatter(tshift(ts)*ones(size(aspk_phase)),aspk_phase,'.k');hold on
            [conva(:,ts),atmap,~] = circle_map(aspk_phase,apath',at',kappa,mapAxis); 
            ainfo(ts) = (atmap/sum(atmap).*conva(:,ts)/(conva(:,ts)'*atmap/sum(atmap)))'*log2(conva(:,ts)/(conva(:,ts)'*atmap/sum(atmap))); %spatial information (bit/spike)
        end
%         pause
        
        [p,pidx] = max(ainfo);
        apeak = [apeak;p];
        apeakt = [apeakt;tshift(pidx)];
        askew = [askew;myskewness(tshift,0,[],ainfo)];       
        
        figure(1)
        imagesc(tshift,mapAxis,conva);axis xy;colorbar;colormap hot
        xlabel('times shift (sec)');ylabel('track position (rads)');
        title(strcat('Convolution of Ambiguous Laps: Cell ',num2str(cells(nc))));
        hold on;plot([0 0],[min(mapAxis) max(mapAxis)],'b');plot([min(tshift) max(tshift)],[0 0],'b');hold off
        figImage = sprintf('%s%s%s%s',cd,strcat('\Convolution_amaps_first','\'),strcat('cell_',num2str(cells(nc))),'.fig');
        bmpImage = sprintf('%s%s%s%s',cd,strcat('\Convolution_amaps_first','\'),strcat('cell_',num2str(cells(nc))),'.bmp');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');   
        plot(tshift,ainfo,'k');
        xlabel('time shift (sec)');
        ylabel('spatial information (bits/spike)');
        title(strcat('Spatial Information Against Spike Shifts: Cell ',num2str(cells(nc))));
        legend('ambiguous laps');
        figImage = sprintf('%s%s%s%s',cd,strcat('\Convolution_ainfo_first','\'),strcat('cell_',num2str(cells(nc))),'first.fig');
        bmpImage = sprintf('%s%s%s%s',cd,strcat('\Convolution_ainfo_first','\'),strcat('cell_',num2str(cells(nc))),'first.bmp');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');             
        close all
    end
end

% data.rpeak = rpeak;
% data.rpeakt = rpeakt;
% data.rskew = rskew;
% 
% filename = sprintf('%s%s%s%s',cd,strcat('\Convolution_rinfo','\'),'data.mat');
% save(filename,'-struct','data');
% 
% data.ppeak = ppeak;
% data.ppeakt = ppeakt;
% data.pskew = pskew;
% 
% filename = sprintf('%s%s%s%s',cd,strcat('\Convolution_pinfo','\'),'data.mat');
% save(filename,'-struct','data');

data.apeak = apeak;
data.apeakt = apeakt;
data.askew = askew;

filename = sprintf('%s%s%s%s',cd,strcat('\Convolution_ainfo_first','\'),'data.mat');
save(filename,'-struct','data');



