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

function Spike_Phase_Hist_MultiU(inFile,s1,s2,f1,f2)

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
    if ii == -1
        ttList = str;
    elseif ii == 0
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

% read the file names from the tt-file list
ttid = fopen(ttList,'r');
jj = 1;
while ~feof(ttid)
    str = fgetl(ttid);
    cells(jj) = {str};
    jj = jj+1;
end
numcells = jj-1;

for ii = 1:numsessions
    disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));
    % Check if subdir for storing images are present. If not, it is
    % created
    dirInfo = dir(sessions{ii});
    found = 0;
    for ss=1:size(dirInfo,1)
        if dirInfo(ss).isdir
            if strcmp(dirInfo(ss).name,strcat('MultiUSpikePhasePlots','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(sessions{ii},strcat('MultiUSpikePhasePlots','\'));
    end
    TS = [];
    % Load data from the .tt and .ncs files, make plots, and store them
    for jj=1:numcells
        disp('Make plots and store them to files');
        disp(sprintf('%s%i',' cell ',jj));
        
        dirFiles = dir(sessions{ii});
        
        if isempty(strmatch(char(cells{jj}),char(dirFiles.name),'exact'))
            continue
        end
        
        tfile = [sessions{ii},cells{jj}];
        [TS] = [TS; loadSpikes(tfile)];
    end
    for jj=1
        filex = [sessions{ii},channels{jj}];
        [samples,ts,tt, Fs, bv, ir] = loadEEG2(filex);
        ch_X = bv*samples;
        filey = [sessions{ii},channels{jj+1}];
        [samples,ts,tt, Fs, bv, ir] = loadEEG2(filey);
        ch_Y = bv*samples;
        %set recording channel against ground and then against average reference
        xcscnum = channels{jj}(4:end-4); % get  tt# corresponding to avgsjj
        xcscnum = str2num(xcscnum);
        xttnum = ceil(xcscnum/4);
        ycscnum = channels{jj+1}(4:end-4); % get  tt# corresponding to avgsjj
        ycscnum = str2num(ycscnum);
        yttnum = ceil(ycscnum/4);
        if ~strcmp(refs{xttnum},'G')
            xrfile = [sessions{ii},refs{xttnum}];
            [samples,ts,tt, Fs, bv, ir] = loadEEG2(xrfile);
            ch_X = ch_X + bv*samples;
        end
        if ~strcmp(refs{yttnum},'G')
            yrfile = [sessions{ii},refs{yttnum}];
            [samples,ts,tt, Fs, bv, ir] = loadEEG2(yrfile);
            ch_Y = ch_Y + bv*samples;
        end
        % slow gamma
        [peak_ind_unique,~,~] = code_for_slow_fast_gamma_indices_windows3_CA1only(ch_X,ch_Y, s1, s2, Fs);
        [~, ~, ~, phaseBin_sPhase2_CA1_cell_ca1_gamma_epochs, phaseBin2_sPhase2_CA1_cell_ca1_gamma_epochs] = CA1_gamma_modulation_of_CA1_spikes(TS, ts, ch_X,peak_ind_unique, s1, s2, Fs);
        if isempty(phaseBin_sPhase2_CA1_cell_ca1_gamma_epochs)
            continue
        end
        data_s.s1 = s1; data_s.s2 = s2; data_s.spikecounts = phaseBin_sPhase2_CA1_cell_ca1_gamma_epochs;   %save data as struct
        filename = sprintf('%s%s%s%s',sessions{ii},strcat('MultiUSpikePhasePlots','\slow_'),cells{jj}(1:end-2),'.mat');
        save(filename,'-struct','data_s');
        figure(1);
        bar([0:30:690],[phaseBin_sPhase2_CA1_cell_ca1_gamma_epochs;phaseBin2_sPhase2_CA1_cell_ca1_gamma_epochs],1)
        xlabel('phase (deg)'); ylabel('spike count')
        title(['Phase Histogram Slow Gamma ',num2str(s1),'-',num2str(s2),' Hz: ',tfile]);
        figImage = sprintf('%s%s%s%s',sessions{ii},strcat('MultiUSpikePhasePlots','\slow_'),cells{jj}(1:end-2),'.fig');
        bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('MultiUSpikePhasePlots','\slow_'),cells{jj}(1:end-2),'.bmp');
        %         epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('SpikePhasePlots','\slow_'),cells{jj}(1:end-2),'.eps');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');
        %         saveas(gcf,epsImage,'eps');
        % fast gamma
        [peak_ind_unique,~,~] = code_for_slow_fast_gamma_indices_windows3_CA1only(ch_X,ch_Y, f1, f2, Fs);
        [~, ~, ~, phaseBin_sPhase2_CA1_cell_ca1_gamma_epochs, phaseBin2_sPhase2_CA1_cell_ca1_gamma_epochs] = CA1_gamma_modulation_of_CA1_spikes(TS, ts, ch_X, peak_ind_unique, f1, f2, Fs);
        if isempty(phaseBin_sPhase2_CA1_cell_ca1_gamma_epochs)
            continue
        end
        data_f.f1 = f1; data_f.f2 = f2; data_f.spikecounts = phaseBin_sPhase2_CA1_cell_ca1_gamma_epochs;   %save data as struct
        filename = sprintf('%s%s%s%s',sessions{ii},strcat('MultiUSpikePhasePlots','\fast_'),cells{jj}(1:end-2),'.mat');
        save(filename,'-struct','data_f');
        figure(2);
        bar(0:30:690,[phaseBin_sPhase2_CA1_cell_ca1_gamma_epochs;phaseBin2_sPhase2_CA1_cell_ca1_gamma_epochs],1)
        xlabel('phase (deg)'); ylabel('spike count')
        title(['Phase Histogram Fast Gamma ',num2str(f1),'-',num2str(f2),' Hz: ',tfile]);
        figImage = sprintf('%s%s%s%s',sessions{ii},strcat('MultiUSpikePhasePlots','\fast_'),cells{jj}(1:end-2),'.fig');
        bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('MultiUSpikePhasePlots','\fast_'),cells{jj}(1:end-2),'.bmp');
        %         epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('SpikePhasePlots','\fast_'),cells{jj}(1:end-2),'.eps');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');
        %         saveas(gcf,epsImage,'eps');
    end
end





