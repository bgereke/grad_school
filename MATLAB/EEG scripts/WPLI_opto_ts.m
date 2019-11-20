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

function WPLI_opto_ts(inFile,freqVec,width,tbins)  % width - number of cycles in wavelet (> 5 advisable); order - whitening

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

%make directory to save figures
savedir = cd;
dirInfo = dir(cd);
found = 0;
for kk=1:size(dirInfo,1)
    if dirInfo(kk).isdir
        if strcmp(dirInfo(kk).name,strcat('WPLIOptoPlots','\'))
            found = 1;
        end
    end
end
if found==0
    mkdir(savedir,strcat('WPLIOptoPlots','\'));
end

for ii = 1%:numchannels
    
    %declare main vars
    outsum_otrig  = zeros(length(freqVec),length(tbins)-1,2);
    outsumW_otrig  = zeros(length(freqVec),length(tbins)-1,2);
    
    disp('Make plots and store them to files');
%     disp(sprintf('%s%i',' CSC ',ii, ' of ',numchannels));
    
    for jj=1:numsessions
        
        %          avref = 0;
        %          for kk=1:numrefs
        %              file = [sessions{jj},avgs{kk}];
        %              [samples,ts,tt, Fs, bv, ir] = loadEEG2(file);
        %              ch_X = bv*samples;
        %              cscnum = avgs{kk}(4:end-4); % get  tt# corresponding to avgsjj
        %              cscnum = str2num(cscnum);
        %              ttnum = ceil(cscnum/4);
        %              if ~strcmp(refs{ttnum},'G')
        %                  rfile = [sessions{jj},refs{ttnum}];
        %                  [samples,ts,tt, Fs, bv, ir] = loadEEG2(rfile);
        %                  ch_X = ch_X + bv*samples;
        %              end
        %              avref = avref + ch_X/numrefs;
        %              clear samples ch_X
        %          end
        
        efile = strcat(sessions{jj},'Events.nev');
        FS = [1 0 1 0 1]; EH = 0; EM = 1;
        [TimeStamps, TTLs, EventStrings] = Nlx2MatEV(efile,FS,EH,EM);
        TTLs(1) = []; TimeStamps(1) = []; TimeStamps = TimeStamps/1000000;
        
        file = [sessions{jj},channels{ii}];
        [samples,ts,tt, Fs, bv, ir] = loadEEG2(file);
        ch_X = bv*samples;
        file = [sessions{jj},channels{ii+1}];
        [samples,ts,tt, Fs, bv, ir] = loadEEG2(file);
        ch_Y = bv*samples;
        
        [ImX] = traces2ImX(ch_X,ch_Y,freqVec,Fs,width);
        
        tt(tt<TimeStamps(1))=[];ImX(tt<TimeStamps(1))=[];
        tt(tt>TimeStamps(end))=[];ImX(tt>TimeStamps(end))=[];
        tt = tt - min(TimeStamps);
        TS = TimeStamps - min(TimeStamps);
        
        for ts=2:length(TS)
            tt_locked = tt-TS(ts-1);
            for t=1:length(tbins)-1
                tidx = tt_locked >= tbins(t) & tt_locked < tbins(t+1);
                if TTLs(ts) == 0
                    outsum_otrig(:,t,1) = outsum_otrig(:,t,1)+sum(ImX(:,tidx),2);
                    outsumW_otrig(:,t,1) = outsumW_otrig(:,t,1)+sum(abs(ImX(:,tidx)),2);
                else
                    outsum_otrig(:,t,2) = outsum_otrig(:,t,2)+sum(ImX(:,tidx),2);
                    outsumW_otrig(:,t,2) = outsumW_otrig(:,t,2)+sum(abs(ImX(:,tidx)),2);
                end
            end
        end
    end
    
    close all
    figure(1);
    
    WPLI_otrig = (outsum_otrig./outsumW_otrig).^2;
    WPLI_plot = [WPLI_otrig(:,:,1), WPLI_otrig(:,:,2), WPLI_otrig(:,:,1)];
    plot_t = [tbins(2:end), tbins(2:end)+max(tbins), tbins(2:end)+2*max(tbins)];
    
    %         imagesc(tt, freqVec, TFR, [0 0.25*max(max(TFR))]); axis xy; xlabel('time (sec)'); ylabel('frequency (Hz)');
    imagesc(plot_t, freqVec, WPLI_plot); axis xy; xlabel('time (sec)'); ylabel('frequency (Hz)');
    hold on; colorbar; colormap hot
    plot([30 30],[min(freqVec) max(freqVec)],'-g');
    plot([60 60],[min(freqVec) max(freqVec)],'-g');
    rectangle('Position',[tbins(1)+max(tbins) max(freqVec)+2,tbins(end),4],'FaceColor','g','EdgeColor','g')
    xlim([tbins(1), 3*tbins(end)]); ylim([min(freqVec) max(freqVec)+4]);
    title(['Opto-triggered WPLI: ',channels{ii}]); colorbar; hold off
    
    data.WPLI_plot = WPLI_plot; data.plot_t = plot_t; data.freqVec = freqVec; data.tbins = tbins;   %save data as struct
    filename = sprintf('%s%s%s%s',savedir,strcat('\WPLIOptoPlots','\'),channels{ii}(1:end-4),channels{ii+1}(1:end-4),'.mat');
    save(filename,'-struct','data');
    figImage = sprintf('%s%s%s%s',savedir,strcat('\WPLIOptoPlots','\'),channels{ii}(1:end-4),channels{ii+1}(1:end-4),'.fig');
    %         bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_'),channels{jj}(1:end-4),'.bmp');
    %         epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_'),channels{jj}(1:end-4),'.eps');
    saveas(gcf,figImage,'fig');
    %         saveas(gcf,bmpImage,'bmp');
    %         saveas(gcf,epsImage,'eps');
    
    
end






