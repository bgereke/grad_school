function bayes_decode(inFile,freqVec,phasebins,width)

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

data.eeg = []; data.TFR = []; data.tsi = []; data.timeVec = []; data.freqVec = []; 
data.bp = []; data.thetadelta = [];

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
    [TFR,timeVec,data.freqVec] = traces2TFR([eeg eeg],freqVec,Fs,width);  
    
    theta = mean(TFR(freqVec>6 & freqVec<12,:),1); %theta
    delta = mean(TFR(freqVec>1 & freqVec<=4,:),1);%delta
    thetadelta = smooth(theta./delta,1000);    
    bp = fftbandpass(eeg,2000,4,6,12,14);%theta
    
    data.eeg = [data.eeg; eeg];
    data.TFR = [data.TFR TFR];
    data.thetadelta = [data.thetadelta;thetadelta];
    data.tsi = [data.tsi;tsi];
    if size(bp,1) == 1
        data.bp = [data.bp bp];
    else
        data.bp = [data.bp;bp];
    end
    data.timeVec = [data.timeVec timeVec];
      
end

%     data.eeg = data.eeg/numchannels;
%     data.TFR = data.TFR/numchannels;

%zscore prior to theta filtering
data.TFRz = zscore(data.TFR,0,2);
data.sg = mean(data.TFR(freqVec>25&freqVec<50,:));
data.sgz = zscore(data.sg);
data.fg = mean(data.TFR(freqVec>55&freqVec<100,:));
data.fgz = zscore(data.fg);

[data] = thetawindows(data);
'still working1'
%for computing the global PFR on data zscored prior to theta filtering
%     for p=1:size(phasebins,2)-1
%         pidx = data.theta_phase >= phasebins(p) & data.theta_phase < phasebins(p+1);
%         data.PFR(:,p) = mean(data.TFRz(:,pidx),2);
%     end

%zscore post theta filtering
tdTFR = [];
for ww = 1:size(data.w,1)
    tdTFR = [tdTFR data.TFR(:,data.widx(ww,1):data.widx(ww,2))];
end
tdsg = mean(tdTFR(freqVec>25&freqVec<50,:),1);
tdfg = mean(tdTFR(freqVec>55&freqVec<100,:),1);
tdsgz = zscore(tdsg,0,2);
tdfgz = zscore(tdfg,0,2);
tdTFRz = zscore(tdTFR,0,2);

%compute tfr's, pfr's, and gamma's for each window
data.wtfr = cell(size(data.w,1),1); %TFR's for each window
data.wtfrz = cell(size(data.w,1),1); %TFRz's for each window
data.wphase = cell(size(data.w,1),1); %theta phase time series for each window
data.wpfr = cell(size(data.w,1),1); %PFR's for each window
data.wsg = cell(size(data.w,1),1); %slow gamma for each window
data.wfg = cell(size(data.w,1),1); %fast gamma for each window
data.wsgz = cell(size(data.w,1),1); %zscored slow gamma for each window
data.wfgz = cell(size(data.w,1),1); %zscored fast gamma for each window
data.PFR = []; %for computing global PFR on theta filtered data
data.PFRv = []; %for computing global PFR variance on theta filtered data
data.wsgmu = nan(size(data.w,1),1); %mean slow gamma for each window
data.wfgmu = nan(size(data.w,1),1); %mean fast gamma for each window

for ww = 1:size(data.w,1)
    data.TFRz(:,data.widx(ww,1):data.widx(ww,2)) = tdTFRz(:,1:data.widx(ww,2)-data.widx(ww,1)+1);
    data.sg(:,data.widx(ww,1):data.widx(ww,2)) = tdsg(:,1:data.widx(ww,2)-data.widx(ww,1)+1);
    data.fg(:,data.widx(ww,1):data.widx(ww,2)) = tdfg(:,1:data.widx(ww,2)-data.widx(ww,1)+1);
    data.sgz(:,data.widx(ww,1):data.widx(ww,2)) = tdsgz(:,1:data.widx(ww,2)-data.widx(ww,1)+1);
    data.fgz(:,data.widx(ww,1):data.widx(ww,2)) = tdfgz(:,1:data.widx(ww,2)-data.widx(ww,1)+1);
    
    data.wsgmu(ww) = mean(tdsg(:,1:data.widx(ww,2)-data.widx(ww,1)+1));
    data.wfgmu(ww) = mean(tdfg(:,1:data.widx(ww,2)-data.widx(ww,1)+1));
    
    tdTFRz(:,1:data.widx(ww,2)-data.widx(ww,1)+1) = [];
    tdsg(:,1:data.widx(ww,2)-data.widx(ww,1)+1) = [];
    tdfg(:,1:data.widx(ww,2)-data.widx(ww,1)+1) = [];
    tdsgz(:,1:data.widx(ww,2)-data.widx(ww,1)+1) = [];
    tdfgz(:,1:data.widx(ww,2)-data.widx(ww,1)+1) = [];
end

data.wsgzmu = zscore(data.wsgmu); %zscored mean slow gamma for each window
data.wfgzmu = zscore(data.wfgmu); %zscored mean fast gamma for each window

for ww = 1:size(data.w,1)
    data.wtfr{ww} = data.TFR(:,data.widx(ww,1):data.widx(ww,2));
    data.wtfrz{ww} = data.TFRz(:,data.widx(ww,1):data.widx(ww,2));
    data.wsg{ww} = data.sg(data.widx(ww,1):data.widx(ww,2));
    data.wsgz{ww} = data.sgz(data.widx(ww,1):data.widx(ww,2));
    data.wfg{ww} = data.fg(data.widx(ww,1):data.widx(ww,2));
    data.wfgz{ww} = data.fgz(data.widx(ww,1):data.widx(ww,2));
    data.wphase{ww} = data.theta_phase(data.widx(ww,1):data.widx(ww,2));
    for p=1:size(phasebins,2)-1
        pidx = data.wphase{ww} >= phasebins(p) & data.wphase{ww} < phasebins(p+1);
        data.wpfr{ww}(:,p) = mean(data.wtfrz{ww}(:,pidx),2);
    end
    data.PFR(:,:,ww) = data.wpfr{ww};
end

data.PFRv = nanvar(data.PFR,[],3);
data.PFR = nanmean(data.PFR,3);
    
'stillworking2'
%decode all the windows
[data] = decodewindows(sessions,ttList,data);   

% Check if subdir for storing images are present. If not, it is
% created
dirInfo = dir(sessions{1}(1:end-7));
found = 0;
for ss=1:size(dirInfo,1)
    if dirInfo(ss).isdir
        if strcmp(dirInfo(ss).name,strcat('DecodingPlots','\'))
            found = 1;
        end
    end
end
if found==0
    mkdir(sessions{1}(1:end-7),strcat('DecodingPlots','\'));
end

%save all data to session folder
filename = sprintf('%s%s%s%s',strcat('DecodingPlots','\'),'cycle_data_thetamins_posttdz_csc1_sepruns_noprior_sessionscombined','.mat');
save(filename,'-struct','data');

clear data
    
    
    
         
        



