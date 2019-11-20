clear
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
        if strcmp(dirInfo(ss).name,strcat('Probability_rmaps','\'))
            found1 = 1;
        end
        if strcmp(dirInfo(ss).name,strcat('Probability_pmaps','\'))
            found2 = 1;
        end     
        if strcmp(dirInfo(ss).name,strcat('Probability_amaps_first','\'))
            found3 = 1;
        end
        if strcmp(dirInfo(ss).name,strcat('Probability_rinfo','\'))
            found4 = 1;
        end    
        if strcmp(dirInfo(ss).name,strcat('Probability_pinfo','\'))
            found5 = 1;
        end
        if strcmp(dirInfo(ss).name,strcat('Probability_ainfo_first','\'))
            found6 = 1;
        end  
    end
end
if found1==0
    mkdir(dd,strcat('\Probability_rmaps','\'));
end
if found2==0
    mkdir(dd,strcat('\Probability_pmaps','\'));
end
if found3==0
    mkdir(dd,strcat('\Probability_amaps_first','\'));
end
if found4==0
    mkdir(dd,strcat('\Probability_rinfo','\'));
end
if found5==0
    mkdir(dd,strcat('\Probability_pinfo','\'));
end
if found6==0
    mkdir(dd,strcat('\Probability_ainfo_first','\'));
end

cells = unique(cmodes(:,3));
numcells = length(cells);
tbins = linspace(-2,2,200);      %time shifts for convolution
D = 100;    %cm
binsize = 5; %cm
% mapAxis = linspace(-pi,pi,pi*D/binsize);
xbins = linspace(-2,2,200);
dx = xbins(2)-xbins(1);
kappa = 400; % Smoothing factor for track location
invh = 10; % Smoothing factor for time shift
pcut = 30;  %number of passes required to compute rate maps
acut = 25;
peakt = [];
peakx = [];
skewt = [];
skewx = [];

for nc = 1:numcells
    
    numr = sum(cmodes(:,3)==cells(nc)&cmodes(:,1)==1);    %number of mode-triggered passes
    nump = sum(cmodes(:,3)==cells(nc)&cmodes(:,1)==-1);
    numa = sum(cmodes(:,3)==cells(nc)&isnan(cmodes(:,1)));    
    
    if numa+numr+nump >= acut 
        apasses = find(cmodes(:,3)==cells(nc));%&isnan(cmodes(:,1)));
        ats = [];   %spike times
        apath = []; %animal trajectories 
%         at = [];    %tracking times
        atshift = []; %tracking times relative to spike time
        
        for ap = apasses'
            aRows = ismember(pass_ts(:,2:6),cmodes(ap,2:6),'rows');
            idx = find(aRows,1,'first'); %first spike only
            ats = [ats;pass_ts(idx,1)];
            aRows = ismember(pass_path(:,2:6),cmodes(ap,2:6),'rows');
            apath = [apath;pass_path(aRows,1)];
%             at = [at;pass_t(aRows,1)];
            atshift = [atshift;pass_t(aRows,1)-pass_ts(idx,1)];
        end  
        
        [PDxt,Mx,Mt] = PDmap(apath,atshift,invh,kappa,xbins,tbins);
        [~,xidx] = max(Mx);
        [~,tidx] = max(Mt);
        peakt = [peakt;tbins(tidx)];
        peakx = [peakx;xbins(xidx)];
        skewt = [skewt;myskewness(tbins,0,[],Mt)];
        skewx = [skewx;myskewness(xbins,0,[],Mx)];
        
        h1 = subplot(2,2,1);        
        imagesc(tbins,xbins,PDxt);axis xy;colormap hot;colorbar;hold on
        scatter(atshift,apath,0.2,'.k');xlim([min(tbins) max(tbins)]);ylim([min(xbins) max(xbins)]);        
        hold on;plot([0 0],[min(xbins) max(xbins)],'b');plot([min(tbins) max(tbins)],[0 0],'b');hold off
        set(h1,'Box','on','YAxisLocation','right');
        xlabel('times shift (sec)');ylabel('track position (rads)');
        title(strcat('Relative Spike Time Probability Distribution: Cell ',num2str(cells(nc)))); 
        h2 = subplot(2,2,2);
        plot(Mx,xbins,'k');xlabel('peak probability')
        set(h2,'Box','on');
        h3 = subplot(2,2,3);
        plot(tbins,Mt,'k');ylabel('spatial information (bits/spike)')  
        set(h3,'YDir','reverse','XAxisLocation','top','Box','on');        
        set(gcf,'color','w');
        figImage = sprintf('%s%s%s%s',cd,strcat('\Probability_amaps_first','\'),strcat('cell_',num2str(cells(nc))),'.fig');
        bmpImage = sprintf('%s%s%s%s',cd,strcat('\Probability_amaps_first','\'),strcat('cell_',num2str(cells(nc))),'.bmp');  
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

data.peakt = peakt;
data.peakx = peakx;
data.skewt = skewt;
data.skewx = skewx;

filename = sprintf('%s%s%s%s',cd,strcat('\Probability_ainfo_first','\'),'data.mat');
save(filename,'-struct','data');



