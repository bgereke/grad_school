% CFC('inFile.txt')
%
% CFC generates cross frequency coherences, plots them and stores them to images and .mat files.
% Multiple sessions will be read from the CSC's specififed in
% 'CSCList.txt'. 'TTList.txt' is necessary for spike data scripts that use
% the same 'infile.txt'.
%
% The input file must be on the following format.
%
% C:\Data\TTList.txt
% C:\Data\CSCList.txt
% C:\Data\Begin 1
% C:\Data\Begin 2
% C:\Data\Begin 3
% C:\Data\Begin 4
% and so on ...
%
% 'CSCList.txt' contains a list of the Neuralynx .csc files to be analyzed.
% All plots will be stored to both bmp and eps imagefiles to a subdirectory in
% the data folder called CFC_plots.

function TFR_w_light(inFile,freqVec,width)  % width - number of cycles in wavelet (> 5 advisable); order - whitening

img_text = 'on';

fid = fopen(inFile,'r');
if fid == -1
    msgbox('Could not open the input file! Make sure the filname and path are correct.','ERROR');
end

% Get sessions and csc-file list from input file
fid = fopen(inFile,'r');
ii = -1;     
numsessions = 0;
while ~feof(fid)
    str = fgetl(fid);
    if ii == 0
        cscList = str;
    elseif ii == 1
        refList = str;
        elseif ii == 2
        ch4avg = str;
    elseif ii > 2
        numsessions  = numsessions+1;
        if ~strcmp(str(end),'\')
            str = strcat(str,'\');
        end
        sessions(numsessions) = {str};        
    end
    ii = ii+1;
end    

% read the file names from the csc-file list
cscid = fopen(cscList,'r');
jj = 1;
while ~feof(cscid)
       str = fgetl(cscid);
       channels(jj) = {str};
       jj = jj+1;
end
numchannels = jj-1;
cscid = fclose('all');

% read the file names of the references for the channels to be used for
% averaging
refid = fopen(refList,'r');
jj = 1;
while ~feof(refid)
       str = fgetl(refid);
       refs(jj) = {str};
       jj = jj+1;
end

% read the file names of the channels to be used for averaging
avgid = fopen(ch4avg,'r');
jj = 1;
while ~feof(avgid)
       str = fgetl(avgid);
       avgs(jj) = {str};
       jj = jj+1;
end
numrefs = jj-1;

for ii = 1:numsessions
    disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));
    % Check if subdir for storing images are present. If not, it is
    % created
    dirInfo = dir(sessions{ii});
    found = 0;
    for kk=1:size(dirInfo,1)
        if dirInfo(kk).isdir
            if strcmp(dirInfo(kk).name,strcat('TFRplots','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(sessions{ii},strcat('TFRplots','\'));
    end
    
    avref = 0;
    for jj=1:numrefs
        file = [sessions{ii},avgs{jj}];
        [samples,ts,tt, Fs, bv, ir] = loadEEG2(file);
        ch_X = bv*samples;
        cscnum = avgs{jj}(4:end-4); % get  tt# corresponding to avgsjj 
        cscnum = str2num(cscnum);   
        ttnum = ceil(cscnum/4);
        if ~strcmp(refs{ttnum},'G')
            rfile = [sessions{ii},refs{ttnum}];
            [samples,ts,tt, Fs, bv, ir] = loadEEG2(rfile);
            ch_X = ch_X + bv*samples;
        end
        avref = avref + ch_X/numrefs;
        clear samples ch_X
    end
    
    efile = strcat(sessions{ii},'Events.nev');
    FS = [1 0 1 0 1]; EH = 0; EM = 1;
    [TimeStamps, TTLs, EventStrings] = Nlx2MatEV(efile,FS,EH,EM);
    TTLs(1) = []; TimeStamps(1) = []; TimeStamps = TimeStamps/1000000;
    
    % Load data from the .ncs files, make plots, and store them
    for jj=1:numchannels
        disp('Make plots and store them to files');
        disp(sprintf('%s%i',' CSC ',jj, ' of ',numchannels));        
        file = [sessions{ii},channels{jj}];
        [samples,ts,tt, Fs, bv, ir] = loadEEG2(file);
        ch_X = bv*samples;
        %set recording channel against ground and then against average reference
        cscnum = channels{jj}(4:end-4); % get  tt# corresponding to avgsjj 
        cscnum = str2num(cscnum);   
        ttnum = ceil(cscnum/4);
        if ~strcmp(refs{ttnum},'G')
            rfile = [sessions{ii},refs{ttnum}];
            [samples,ts,tt, Fs, bv, ir] = loadEEG2(rfile);
            ch_X = ch_X + bv*samples;            
        end
        ch_X = ch_X - avref;  
        clear samples
        TS = TimeStamps - min(tt);
        tt = tt - min(tt);
%         rfile = [sessions{ii},'R2.ncs'];
%         [samples,ts,tt, Fs, bv, ir] = loadEEG2(rfile);
%         Ref = bv*samples;
%         ch_X = ch_X + Ref;
%         ch_Xw = prewhitening(ch_X,3);   % prewhiten signal for better visualisation
        [TFR,timeVec,freqVec] = traces2TFR([ch_X ch_X],freqVec,Fs,width);
        
%         TFRperc = zeros(size(TFR));
% 
%         for i=1:size(TFR,1)
%             TFRperc(i,:) = tiedrank(TFR(i,:))/size(TFR,2);
%         end
%         theta = mean(TFR(freqVec>6 & freqVec<12,:),1); %theta
%         delta = mean(TFR(freqVec>1 & freqVec<=4,:),1);%delta
%         thetadelta = smooth(theta./delta,1000);
        TFRz = zscore(TFR,0,2);
        data.tt = tt; data.freqvec = freqVec; data.tfr = TFR; %data.td = thetadelta;  %save data as struct
        filename = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\'),channels{jj}(1:end-4),'.mat');
        save(filename,'-struct','data');
        close all; figure(1); hold on        
%         imagesc(tt, freqVec, TFR, [0 0.25*max(max(TFR))]); axis xy; xlabel('time (sec)'); ylabel('frequency (Hz)');
        imagesc(tt, freqVec, TFRz, [-2 2]); axis xy; xlabel('time (sec)'); ylabel('frequency (Hz)');
        
        if TTLs(1) == 1            
            rectangle('Position',[min(tt) max(freqVec)+2,TS(1),4],'FaceColor','g','EdgeColor','g')
        end
        for i = 2:length(TTLs)            
            if TTLs(i) == 1  
                rectangle('Position',[TS(i-1),max(freqVec)+2,TS(i)-TS(i-1),4],'FaceColor','g','EdgeColor','g')
            end
            tidx = length(tt(tt<=TS(i)))+1;
        end        
        if TTLs(end) == 0
            rectangle('Position',[TS(end),max(freqVec)+2,max(tt)-TS(end),4],'FaceColor','g','EdgeColor','g')
        end
        xlim([min(tt) max(tt)]); ylim([min(freqVec) max(freqVec)+4]);
        title(['Time Frequency Representation: ',file]); colorbar;
        figImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_'),channels{jj}(1:end-4),'.fig');
%         bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_'),channels{jj}(1:end-4),'.bmp');
%         epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_'),channels{jj}(1:end-4),'.eps');
        saveas(gcf,figImage,'fig');
%         saveas(gcf,bmpImage,'bmp');
%         saveas(gcf,epsImage,'eps');
        clear TFR TFRz data
    end
end




