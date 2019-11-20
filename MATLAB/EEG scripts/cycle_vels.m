function cycle_vels(inFile,freqVec,phasebins,width)

img_text = 'on';

fid = fopen(inFile,'r');
if fid == -1
    msgbox('Could not open the input file! Make sure the filname and path are correct.','ERROR');
end

% Get sessions and csc-file list from input file
fid = fopen(inFile,'r');
ii = -1;     
while ~feof(fid)
    str = fgetl(fid);
    if ii == -1
        ttList = str;
    elseif ii == 0
        cscList = str;
    elseif ii > 0
        if ~strcmp(str(end),'\')
            str = strcat(str,'\');
        end
        sessions(ii) = {str};
    end
    ii = ii+1;
end
numsessions = ii-1;   

% read the file names from the tt-file list
ttid = fopen(ttList,'r');
jj = 1;
while ~feof(ttid)
       str = fgetl(ttid);
       cells(jj) = {str};
       jj = jj+1;
end
numcells = jj-1;

% read the file names from the csc-file list
cscid = fopen(cscList,'r');
kk = 1;
while ~feof(cscid)
       str = fgetl(cscid);
       channels(kk) = {str};
       kk = kk+1;
end
numchannels = kk-1;

% data.eeg = []; data.TFR = []; data.tsi = [];% data.timeVec = []; data.freqVec = []; 
% tbounds = []; %[min(beg1) max(beg1);min(beg2) max(beg2);...]
%data.bp = []; data.thetadelta = [];

f1 = 6; f2 = 12; %theta limits
f3 = 1; f4 = 4; %delta limits
data.wtfrz = []; data.w = []; data.tsi = []; data.wphase = []; data.wthetadelta = [];

for ii = 1:numsessions
    disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));   
    
    %compute global TFR, PFR, and get theta windows
%     data.eeg = 0;
%     data.TFR = 0;
    %for averaging the eeg across tetrodes
%     for kk = 1:numchannels
%         file = [sessions{ii},channels{kk}];
%         [ch_X,tsi,ts,Fs] = loadEeg8(file);
%         data.eeg = data.eeg + ch_X;
%         [tfr,data.timeVec,data.freqVec] = traces2TFR([ch_X ch_X],freqVec,Fs,width);
%         data.TFR = data.TFR + tfr;
%         data.tsi = tsi;
%     end

    %for specifying a specific csc to use
    file = [sessions{ii},channels{1}];
    [eeg,tsi,ts,Fs] = loadEeg8(file); 
%     rfile = [sessions{ii},'R2.ncs'];
%     [samples,rts,rtt, rFs, bv, ir] = loadEEG2(rfile);
%     Ref = bv*samples;
%     eeg = eeg + Ref;
    
    [TFR,timeVec,data.freqVec] = traces2TFR([eeg eeg],freqVec,Fs,width);
    
    bp = fftbandpass(eeg,2000,f1-2,f1,f2,f2+2);%theta
    DTAS = hilbert(bp);
    phase = angle(DTAS);
    
    thetadelta = mean(TFR(data.freqVec>=f1&data.freqVec<=f2,:))./...
        mean(TFR(data.freqVec>=f3&data.freqVec<=f4,:));
    
    TFRz = zscore(TFR,0,2); clear TFR
    size_w = 0.2;
    w = [(tsi(1):size_w:tsi(end)-size_w)' (tsi(1)+size_w:size_w:tsi(end))'];
    wtfrz = cell(1,size(w,1));
    wphase = cell(1,size(w,1));
    wthetadelta = cell(1,size(w,1));
    widx = zeros(size(w));
    widx(1,1) = 1;
    for i=1:size(w,1)
        [~,widx(i,2)] = min((tsi - w(i,2)).^2);
        widx(i+1,1) = widx(i,2);
        
        %wtfr{i} = TFR(:,widx(i,1):widx(i,2));
        wtfrz{i} = TFRz(:,widx(i,1):widx(i,2));
        wphase{i} = phase(widx(i,1):widx(i,2));
        wthetadelta{i} = thetadelta(widx(i,1):widx(i,2));
    end
    
    data.w = [data.w; w]; data.wtfrz = [data.wtfrz wtfrz]; data.tsi = [data.tsi;tsi];
    data.wphase = [data.wphase wphase]; data.wthetadelta = [data.wthetadelta wthetadelta];
    clear TFRz w wtfrz widx tsi eeg bp phase wphase thetadelta wthetadelta
%     
%     theta = mean(TFR(freqVec>6 & freqVec<12,:),1); %theta
%     delta = mean(TFR(freqVec>1 & freqVec<=4,:),1);%delta
%     thetadelta = smooth(theta./delta,1000);    
%     bp = fftbandpass(eeg,2000,4,6,12,14);%theta
    
%     data.eeg = [data.eeg; eeg];
%     data.TFR = [data.TFR TFR];
% %     data.thetadelta = [data.thetadelta;thetadelta];
%     data.tsi = [data.tsi;tsi];
%     tbounds = [tbounds;tsi(1) tsi(end)];
%     if size(bp,1) == 1
%         data.bp = [data.bp bp];
%     else
%         data.bp = [data.bp;bp];
%     end
%     data.timeVec = [data.timeVec timeVec];
      
end

%     data.eeg = data.eeg/numchannels;
%     data.TFR = data.TFR/numchannels;

%zscore prior to theta filtering
% data.TFRz = zscore(data.TFR,0,2);
% 
% %[data] = thetawindows(data);
% 
% size_w = 0.2; %temporal bin over which to compute velocity
% data.w = nan(ceil((max(data.tsi)-min(data.tsi))/size_w),2);
% data.widx = data.w;
% data.w(1,1) = data.tsi(1); 
% data.widx(1,1) = 1; 
% for i=1:size(data.w,1)-1
%     
%     if data.w(i,1) > tbounds(1,1) && (data.w(i,1)+size_w) <= tbounds(1,2)
%         [~,data.widx(i,2)] = min((data.tsi - (data.w(i,1)+size_w)).^2);
%         data.w(i,2) = data.tsi(data.widx(i,2));
%         data.w(i+1,1) = data.w(i,2);
%         data.widx(i+1,1) = data.widx(i,2);
%         if data.widx(i,2) + size_w > max(data.tsi)
%             break
%         end
%     elseif data.w(i,1) > tbounds(2,1) && (data.w(i,1)+size_w) <= tbounds(2,2)
%         [~,data.widx(i,2)] = min((data.tsi - (data.w(i,1)+size_w)).^2);
%         data.w(i,2) = data.tsi(data.widx(i,2));
%         data.w(i+1,1) = data.w(i,2);
%         data.widx(i+1,1) = data.widx(i,2);
%         if data.widx(i,2) + size_w > max(data.tsi)
%             break
%         end
%     elseif data.w(i,1) > tbounds(3,1) && (data.w(i,1)+size_w) <= tbounds(3,2)
%         [~,data.widx(i,2)] = min((data.tsi - (data.w(i,1)+size_w)).^2);
%         data.w(i,2) = data.tsi(data.widx(i,2));
%         data.w(i+1,1) = data.w(i,2);
%         data.widx(i+1,1) = data.widx(i,2);
%         if data.widx(i,2) + size_w > max(data.tsi)
%             break
%         end
%     end
% end
% data.w = data.w(~any(isnan(data.w),2),:);
% data.widx = data.widx(~any(isnan(data.widx),2),:);
% 'still working1'
%for computing the global PFR on data zscored prior to theta filtering
%     for p=1:size(phasebins,2)-1
%         pidx = data.theta_phase >= phasebins(p) & data.theta_phase < phasebins(p+1);
%         data.PFR(:,p) = mean(data.TFRz(:,pidx),2);
%     end

%zscore post theta filtering
% tdTFR = [];
% for ww = 1:size(data.w,1)
%     tdTFR = [tdTFR data.TFR(:,data.widx(ww,1):data.widx(ww,2))];
% end
% 
% % tdTFRz = (tdTFR - repmat(mean(tdTFR,2),1,size(tdTFR,2)))./repat(std(tdTFR,2),1,size(tdTFR,2));
% 
% %compute tfr's, pfr's, and gamma's for each window
% data.wtfr = cell(size(data.w,1),1); %TFR's for each window
% % data.wtfrz = cell(size(data.w,1),1); %TFRz's for each window
% data.wphase = cell(size(data.w,1),1); %theta phase time series for each window
% % data.wpfr = cell(size(data.w,1),1); %PFR's for each window
% % data.PFR = []; %for computing global PFR on theta filtered data
% % data.PFRv = []; %for computing global PFR variance on theta filtered data
% 
% % for ww = 1:size(data.w,1)
% %     data.TFRz(:,data.widx(ww,1):data.widx(ww,2)) = tdTFRz(:,1:data.widx(ww,2)-data.widx(ww,1)+1);    
% %     tdTFRz(:,1:data.widx(ww,2)-data.widx(ww,1)+1) = [];  
% % end
% 
% 
% for ww = 1:size(data.w,1)
%     data.wtfr{ww} = data.TFR(:,data.widx(ww,1):data.widx(ww,2));
%     data.wtfrz{ww} = data.TFRz(:,data.widx(ww,1):data.widx(ww,2));
% %     data.wphase{ww} = data.theta_phase(data.widx(ww,1):data.widx(ww,2));
% %     for p=1:size(phasebins,2)-1
% %         pidx = data.wphase{ww} >= phasebins(p) & data.wphase{ww} < phasebins(p+1);
% %         data.wpfr{ww}(:,p) = mean(data.wtfrz{ww}(:,pidx),2);
% %     end
% %     data.PFR(:,:,ww) = data.wpfr{ww};
% end

% data.PFRv = nanvar(data.PFR,[],3);
% data.PFR = nanmean(data.PFR,3);
    
'stillworking2'
%decode all the windows
[data] = getvel(sessions,data);   

% Check if subdir for storing images are present. If not, it is
% created
dirInfo = dir(sessions{1}(1:end-7));
found = 0;
for ss=1:size(dirInfo,1)
    if dirInfo(ss).isdir
        if strcmp(dirInfo(ss).name,strcat('VelData','\'))
            found = 1;
        end
    end
end
if found==0
    mkdir(sessions{1}(1:end-7),strcat('VelData','\'));
end
data.w = data.w-min(data.tsi);
data.tsi = data.tsi-min(data.tsi);
%save all data to session folder
filename = sprintf('%s%s%s%s',strcat('VelData','\'),channels{1},'.mat');
save(filename,'-struct','data');

clear data
    
    
    
         
        



