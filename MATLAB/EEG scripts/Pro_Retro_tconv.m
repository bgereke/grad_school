%Spatial Convolution of mode-triggered fields
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
        if strcmp(dirInfo(ss).name,strcat('tConvolution_rmaps','\'))
            found1 = 1;
        end
        if strcmp(dirInfo(ss).name,strcat('tConvolution_pmaps','\'))
            found2 = 1;
        end     
        if strcmp(dirInfo(ss).name,strcat('tConvolution_amaps_last','\'))
            found3 = 1;
        end
        if strcmp(dirInfo(ss).name,strcat('tConvolution_rinfo','\'))
            found4 = 1;
        end    
        if strcmp(dirInfo(ss).name,strcat('tConvolution_pinfo','\'))
            found5 = 1;
        end
        if strcmp(dirInfo(ss).name,strcat('tConvolution_ainfo_last','\'))
            found6 = 1;
        end  
    end
end
if found1==0
    mkdir(dd,strcat('\tConvolution_rmaps','\'));
end
if found2==0
    mkdir(dd,strcat('\tConvolution_pmaps','\'));
end
if found3==0
    mkdir(dd,strcat('\tConvolution_amaps_last','\'));
end
if found4==0
    mkdir(dd,strcat('\tConvolution_rinfo','\'));
end
if found5==0
    mkdir(dd,strcat('\tConvolution_pinfo','\'));
end
if found6==0
    mkdir(dd,strcat('\tConvolution_ainfo_last','\'));
end

cells = unique(cmodes(:,3));
numcells = length(cells);
D = 100;    %cm
binsize = 5; %cm
% mapAxis = linspace(-pi,pi,pi*D/binsize);
tshift = -10:0.02:10;
mapAxis = linspace(-2,2,150);
dx = mapAxis(2)-mapAxis(1);
mapAxis = mapAxis + dx/2; mapAxis(end) = [];
invh = 4; % Smoothing factor when calculating the ratemap
pcut = 100;  %number of passes required to compute rate maps
acut = 25;
rpeaki = [];
rpeakl = [];
rpeakt = [];
rskew = [];
rposi = [];
rnegi = [];
rposd = [];
rnegd = [];
ppeaki = [];
ppeakl = [];
ppeakt = [];
pskew = [];
pposi = [];
pnegi = [];
pposd = [];
pnegd = [];
apeaki = [];
apeakl = [];
apeakt = [];
askew = [];
aposi = [];
anegi = [];
aposd = [];
anegd = [];

for nc = 1:numcells
    
    convr = zeros(length(tshift),length(mapAxis));  %mode-triggered convolution maps
    convp = convr; conva = convr; 
    rinfo = zeros(length(mapAxis),1);     %mode triggered spatial information
    pinfo = rinfo; ainfo = rinfo;
    
    numr = sum(cmodes(:,3)==cells(nc)&cmodes(:,1)==1);  %number of mode-triggered passes
    nump = sum(cmodes(:,3)==cells(nc)&cmodes(:,1)==-1);
    numa = sum(cmodes(:,3)==cells(nc)&isnan(cmodes(:,1)));
    
    if numr >= pcut 
        
        rpasses = find(cmodes(:,3)==cells(nc)&cmodes(:,1)==1);
        rts = cell(length(mapAxis),1);   %mode triggered spike times
        rpath = cell(length(mapAxis),1); %mode triggered animal trajectories 
        RT = cell(length(mapAxis),1);    %mode-triggered tracking times
        rtpath = cell(length(mapAxis),1);       %time from each point on trajectory to nearest track location
        rspk_t = cell(length(mapAxis),1);      %time from each spike to nearest track location
        
        %compute time of spike and trajectory relative to each point on the map
        
        for rp = rpasses'
            %get pass spike times and trajectory
            rRows = ismember(pass_ts(:,2:6),cmodes(rp,2:6),'rows');
            idx = find(rRows,1,'first'); %first spike only
            rts = pass_ts(idx,1);
            rRows = ismember(pass_path(:,2:6),cmodes(rp,2:6),'rows');
            rpath = pass_path(rRows,1);
            rt = pass_t(rRows,1);
            
            for mp = 1:length(mapAxis)
                mpt = at(abs(apath-mapAxis(mp))<=dx/2); %times animal is at map point mp
                if ~isempty(mpt)    %if the animal doesn't visit the location on that pass
                    [~,idx] = min((mpt-rts).^2);
                    rspk_t{mp} = [rspk_t{mp};mpt(idx)-rts];
                    rtpath{mp} = [rtpath{mp};mpt(idx)-rt];
                    RT{mp} = [RT{mp};rt];
                else
                    continue
                end
            end  
        end
%         aspk_t{1,mapAxis>1.5&mapAxis<2}
        for mp = 1:length(mapAxis)
            %compute the relative rate maps and temporal information for each point on the map
            [convr(:,mp),rtmap,~] = tshift_map(rspk_t{mp},rtpath{mp}',RT{mp}',invh,tshift);
%             rinfo(mp) = (rtmap/sum(rtmap).*convr(:,mp)/(convr(:,mp)'*rtmap/sum(rtmap)))'*log2(convr(:,mp)/(convr(:,mp)'*rtmap/sum(rtmap))); %spatial information (bit/spike)
            rinfo(mp) = max(convr(:,mp));
        end
        [p,pidx] = max(rinfo);
        rpeaki = [rpeaki;p];
        rpeakl = [rpeakl;mapAxis(pidx)];
        [~,idx] = max(convr(:,pidx));
        rpeakt = [rpeakt;tshift(idx)];
        [~,negidx] = min(rinfo(1:pidx));
        [~,posidx] = min(rinfo(pidx:end));
        rskew = [rskew;myskewness(mapAxis(negidx:posidx+pidx-1),0,[],rinfo(negidx:posidx+pidx-1))];
        rposi = [rposi;sum(rinfo(pidx:posidx+pidx-1))];
        rnegi = [rnegi;sum(rinfo(negidx:pidx))];
        rposd = [rposd;mapAxis(posidx+pidx-1)-mapAxis(pidx)];
        rnegd = [rnegd;mapAxis(negidx)-mapAxis(pidx)];
        
        figure(1)
%         for mp = 1:length(mapAxis)-1
%             col = linspace(1,10,length(aspk_t{mp}));
%             scatter(mapAxis(mp)*ones(1,length(aspk_t{mp})),aspk_t{mp},[],col)
%             hold on
%         end
%         colormap jet
%         xlabel('track location (rads)');ylabel('time to location (sec)');pause
        plot(mapAxis(1:end),rinfo,'r');hold on;
        xlabel('track location (rads)');
        ylabel('temporal information (bits/spike)');
        title(strcat('Time to Location Information Retrospective Laps: Cell ',num2str(cells(nc))));
        legend('retrospective laps');hold off;
        figImage = sprintf('%s%s%s%s',cd,strcat('\tConvolution_rinfo','\'),strcat('cell_',num2str(cells(nc))),'first.fig');
        bmpImage = sprintf('%s%s%s%s',cd,strcat('\tConvolution_rinfo','\'),strcat('cell_',num2str(cells(nc))),'first.bmp');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');
        imagesc(mapAxis(1:end),tshift,convr);axis xy;colorbar;colormap hot
        xlabel('track location (rads)');ylabel('time to location (sec)');
        title(strcat('Time to Location Map Retrospective Laps: Cell ',num2str(cells(nc))));
        hold on;plot([0 0],[min(tshift) max(tshift)],'b');plot([min(mapAxis) max(mapAxis)],[0 0],'b');hold off
        figImage = sprintf('%s%s%s%s',cd,strcat('\tConvolution_rmaps','\'),strcat('cell_',num2str(cells(nc))),'first.fig');
        bmpImage = sprintf('%s%s%s%s',cd,strcat('\tConvolution_rmaps','\'),strcat('cell_',num2str(cells(nc))),'first.bmp');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');         
        close all
    end
    if nump >= pcut 
        
        ppasses = find(cmodes(:,3)==cells(nc)&cmodes(:,1)==-1);
        pts = cell(length(mapAxis),1);   %mode triggered spike times
        ppath = cell(length(mapAxis),1); %mode triggered animal trajectories 
        PT = cell(length(mapAxis),1);    %mode-triggered tracking times
        ptpath = cell(length(mapAxis),1);       %time from each point on trajectory to nearest track location
        pspk_t = cell(length(mapAxis),1);      %time from each spike to nearest track location
        
        %compute time of spike and trajectory relative to each point on the map
        
        for pp = ppasses'
            %get pass spike times and trajectory
            pRows = ismember(pass_ts(:,2:6),cmodes(pp,2:6),'rows');
            idx = find(pRows,1,'first'); %first spike only
            pts = pass_ts(idx,1);
            pRows = ismember(pass_path(:,2:6),cmodes(pp,2:6),'rows');
            ppath = pass_path(pRows,1);
            pt = pass_t(pRows,1);
            
            for mp = 1:length(mapAxis)
                mpt = at(abs(apath-mapAxis(mp))<=dx/2); %times animal is at map point mp
                if ~isempty(mpt)    %if the animal doesn't visit the location on that pass
                    [~,idx] = min((mpt-pts).^2);
                    pspk_t{mp} = [pspk_t{mp};mpt(idx)-pts];
                    ptpath{mp} = [ptpath{mp};mpt(idx)-pt];
                    PT{mp} = [PT{mp};pt];
                else
                    continue
                end
            end  
        end
%         aspk_t{1,mapAxis>1.5&mapAxis<2}
        for mp = 1:length(mapAxis)
            %compute the relative rate maps and temporal information for each point on the map
            [convp(:,mp),ptmap,~] = tshift_map(pspk_t{mp},ptpath{mp}',PT{mp}',invh,tshift);
%             pinfo(mp) = (ptmap/sum(ptmap).*convp(:,mp)/(convp(:,mp)'*ptmap/sum(ptmap)))'*log2(convp(:,mp)/(convp(:,mp)'*ptmap/sum(ptmap))); %spatial information (bit/spike)
            pinfo(mp) = max(convp(:,mp));
        end
        [p,pidx] = max(pinfo);
        ppeaki = [ppeaki;p];
        ppeakl = [ppeakl;mapAxis(pidx)];
        [~,idx] = max(convp(:,pidx));
        ppeakt = [ppeakt;tshift(idx)];
        [~,negidx] = min(pinfo(1:pidx));
        [~,posidx] = min(pinfo(pidx:end));
        pskew = [pskew;myskewness(mapAxis(negidx:posidx+pidx-1),0,[],pinfo(negidx:posidx+pidx-1))]; 
        pposi = [pposi;sum(pinfo(pidx:posidx+pidx-1))];
        pnegi = [pnegi;sum(pinfo(negidx:pidx))];
        pposd = [pposd;mapAxis(posidx+pidx-1)-mapAxis(pidx)];
        pnegd = [pnegd;mapAxis(negidx)-mapAxis(pidx)];
        
        figure(1)
%         for mp = 1:length(mapAxis)-1
%             col = linspace(1,10,length(aspk_t{mp}));
%             scatter(mapAxis(mp)*ones(1,length(aspk_t{mp})),aspk_t{mp},[],col)
%             hold on
%         end
%         colormap jet
%         xlabel('track location (rads)');ylabel('time to location (sec)');pause
        plot(mapAxis(1:end),pinfo,'b');hold on;
        xlabel('track location (rads)');
        ylabel('temporal information (bits/spike)');
        title(strcat('Time to Location Information Prospective Laps: Cell ',num2str(cells(nc))));
        legend('prospective laps');hold off;
        figImage = sprintf('%s%s%s%s',cd,strcat('\tConvolution_pinfo','\'),strcat('cell_',num2str(cells(nc))),'first.fig');
        bmpImage = sprintf('%s%s%s%s',cd,strcat('\tConvolution_pinfo','\'),strcat('cell_',num2str(cells(nc))),'first.bmp');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');
        imagesc(mapAxis(1:end),tshift,convp);axis xy;colorbar;colormap hot
        xlabel('track location (rads)');ylabel('time to location (sec)');
        title(strcat('Time to Location Map Prospective Laps: Cell ',num2str(cells(nc))));
        hold on;plot([0 0],[min(tshift) max(tshift)],'b');plot([min(mapAxis) max(mapAxis)],[0 0],'b');hold off
        figImage = sprintf('%s%s%s%s',cd,strcat('\tConvolution_pmaps','\'),strcat('cell_',num2str(cells(nc))),'first.fig');
        bmpImage = sprintf('%s%s%s%s',cd,strcat('\tConvolution_pmaps','\'),strcat('cell_',num2str(cells(nc))),'first.bmp');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');        
        close all
    end
    if numa+numr+nump >= acut 
        apasses = find(cmodes(:,3)==cells(nc));%&isnan(cmodes(:,1)));
        ats = cell(length(mapAxis),1);   %mode triggered spike times
        apath = cell(length(mapAxis),1); %mode triggered animal trajectories 
        AT = cell(length(mapAxis),1);    %mode-triggered tracking times
        atpath = cell(length(mapAxis),1);       %time from each point on trajectory to nearest track location
        aspk_t = cell(length(mapAxis),1);      %time from each spike to nearest track location
        
        %compute time of spike and trajectory relative to each point on the map
        
        for ap = apasses'
            %get pass spike times and trajectory
            aRows = ismember(pass_ts(:,2:6),cmodes(ap,2:6),'rows');
            idx = find(aRows,1,'last'); %first spike only
            ats = pass_ts(idx,1);
            aRows = ismember(pass_path(:,2:6),cmodes(ap,2:6),'rows');
            apath = pass_path(aRows,1);
            at = pass_t(aRows,1);
            
            for mp = 1:length(mapAxis)
                mpt = at(abs(apath-mapAxis(mp))<=dx/2); %times animal is at map point mp
                if ~isempty(mpt)    %if the animal doesn't visit the location on that pass
                    [~,idx] = min((mpt-ats).^2);
                    aspk_t{mp} = [aspk_t{mp};mpt(idx)-ats];                    
%                     atpath{mp} = [atpath{mp};mpt(idx)-at];                    
%                     AT{mp} = [AT{mp};at];
                else
                    continue
                end
            end  
        end
%         figure(2)
%         plot(at,apath);ylim([-2 2])
%         aspk_t{1,mapAxis>1.5&mapAxis<2}
        for mp = 1:length(mapAxis)
%             %compute the relative rate maps and temporal information for each point on the map
            [conva(:,mp),~] = tshift_map(aspk_t{mp},invh,tshift);
%             ainfo(mp) = (1/length(tshift)*conva(:,mp)/sum(conva(:,mp)'/length(tshift)))'*log2(conva(:,mp)/sum(conva(:,mp)'/length(tshift))); %temporal information (bit/spike)
            ainfo(mp) = 1/std(aspk_t{mp});
        end
        conva = conva/length(apasses);
%         ainfo = max(conva);
        [p,pidx] = max(ainfo);
        apeaki = [apeaki;p];
        apeakl = [apeakl;mapAxis(pidx)];
        [~,idx] = max(conva(:,pidx));
        apeakt = [apeakt;tshift(idx)];
        [~,negidx] = min(ainfo(1:pidx));
        [~,posidx] = min(ainfo(pidx:end));
        askew = [askew;myskewness(mapAxis(negidx:posidx+pidx-1),0,[],ainfo(negidx:posidx+pidx-1))];  
        aposi = [aposi;sum(ainfo(pidx:posidx+pidx-1))];
        anegi = [anegi;sum(ainfo(negidx:pidx))];
        aposd = [aposd;mapAxis(posidx+pidx-1)-mapAxis(pidx)];
        anegd = [anegd;mapAxis(negidx)-mapAxis(pidx)];
        
        figure(1)
%         for mp = 1:length(mapAxis)-1
%             col = linspace(1,10,length(aspk_t{mp}));
% %             scatter(mapAxis(mp)*ones(1,length(aspk_t{mp})),aspk_t{mp},[],col)
%             scatter(mapAxis(mp)*ones(1,length(aspk_t{mp})),aspk_t{mp},'.k')
%             hold on
%         end
%         colormap jet
%         hold off
%         xlabel('track location (rads)');ylabel('time to location (sec)');pause
        plot(mapAxis(1:end),ainfo,'k');hold on;
        xlabel('track location (rads)');
        ylabel('temporal information (bits/spike)');
        title(strcat('Time to Location Information Ambiguous Laps: Cell ',num2str(cells(nc))));
        legend('ambiguous laps');hold off;
        figImage = sprintf('%s%s%s%s',cd,strcat('\tConvolution_ainfo_last','\'),strcat('cell_',num2str(cells(nc))),'varlast.fig');
        bmpImage = sprintf('%s%s%s%s',cd,strcat('\tConvolution_ainfo_last','\'),strcat('cell_',num2str(cells(nc))),'varlast.bmp');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');
%         imagesc(mapAxis(1:end),tshift,conva);axis xy;colorbar;colormap hot
%         xlabel('track location (rads)');ylabel('time to location (sec)');
%         title(strcat('Time to Location Map Ambiguous Laps: Cell ',num2str(cells(nc))));
%         hold on;plot([0 0],[min(tshift) max(tshift)],'b');plot([min(mapAxis) max(mapAxis)],[0 0],'b');hold off        
%         figImage = sprintf('%s%s%s%s',cd,strcat('\tConvolution_amaps_last','\'),strcat('cell_',num2str(cells(nc))),'last.fig');
%         bmpImage = sprintf('%s%s%s%s',cd,strcat('\tConvolution_amaps_last','\'),strcat('cell_',num2str(cells(nc))),'last.bmp');
%         saveas(gcf,figImage,'fig');
%         saveas(gcf,bmpImage,'bmp');        
        close all
    end
end

% data.rpeaki = rpeaki;
% data.rpeakl = rpeakl;
% data.rpeakt = rpeakt;
% data.rskew = rskew;
% data.rposi = rposi;
% data.rnegi = rnegi;
% data.rposd = rposd;
% data.rnegd = rnegd;
% 
% filename = sprintf('%s%s%s%s',cd,strcat('\tConvolution_rinfo','\first_'),'data.mat');
% save(filename,'-struct','data');
% clear data
% 
% data.ppeaki = ppeaki;
% data.ppeakl = ppeakl;
% data.ppeakt = ppeakt;
% data.pskew = pskew;
% data.pposi = pposi;
% data.pnegi = pnegi;
% data.pposd = pposd;
% data.pnegd = pnegd;
% 
% filename = sprintf('%s%s%s%s',cd,strcat('\tConvolution_pinfo','\first_'),'data.mat');
% save(filename,'-struct','data');
% clear data

data.apeaki = apeaki;
data.apeakl = apeakl;
data.apeakt = apeakt;
data.askew = askew;
data.aposi = aposi;
data.anegi = anegi;
data.aposd = aposd;
data.anegd = anegd;

filename = sprintf('%s%s%s%s',cd,strcat('\tConvolution_ainfo_last','\last_'),'vardata.mat');
save(filename,'-struct','data');
clear data


