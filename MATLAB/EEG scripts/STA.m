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

function STA(inFile)

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

for ii = 1:numsessions
    disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));
    % Check if subdir for storing images are present. If not, it is
    % created
    dirInfo = dir(sessions{ii});
    found = 0;
    for ss=1:size(dirInfo,1)
        if dirInfo(ss).isdir
            if strcmp(dirInfo(ss).name,strcat('STAPlots','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(sessions{ii},strcat('STAPlots','\'));
    end
    
    % Load data from the .tt and .ncs files, make plots, and store them
    for jj=1:numcells
        disp('Make plots and store them to files');
        disp(sprintf('%s%i',' cell ',jj));
            
        dirFiles = dir(sessions{ii});
       
        if isempty(strmatch(char(cells{jj}),char(dirFiles.name),'exact'))
            continue
        end
         
        tfile = [sessions{ii},cells{jj}];
        [TS] = loadSpikes(tfile);            
        
        tnumber = cells{jj}(3); % gets the cell's tetrode
            if strcmp(cells{jj}(4),'_') == 0   % the tetrode number could be 2 digits
                tnumber = strcat(tnumber,cells{jj}(4));
            end
        tnumber = str2num(tnumber);  
        
        cfile = [];
        for kk=1:numchannels
            cscnumber = channels{kk}(4); % makes list of channel numbers
                if strcmp(channels{kk}(5),'.') == 0
                    cscnumber = strcat(cscnumber,channels{kk}(5));
                end
            cscnumber = str2num(cscnumber);
            if cscnumber > 4*(tnumber-1) && cscnumber <= 4*tnumber
                cfile = [sessions{ii},channels{kk}]; break
            end
        end
        
        if isempty(cfile)
            msgbox('Could not find tetrode channel! Make sure CSCList.txt contains a channel corresponding to each of the tetrodes in TTList.txt.','ERROR');
        end
              
        [samples,ts,tt, Fs, bv, ir] = loadEEG2(cfile);
        ch_X = bv*samples;
        
        figure(1);
        [eeg_epochs_grouped, eeg_epoch_index2] = find_spike_times_rs15_STA(TS, ts, ch_X, Fs);
        data.epochs_grouped = eeg_epochs_grouped; data.epoch_idx = eeg_epoch_index2;    %save data as struct
        filename = sprintf('%s%s%s%s',sessions{ii},strcat('STAPlots','\'),cells{jj}(1:end-2),'.mat');
        save(filename,'-struct','data');
        plot((1:length(mean(eeg_epochs_grouped,2)))/Fs*1000,mean(eeg_epochs_grouped,2));
        xlabel('time (ms)'); ylabel('uV')
        title(['Spike Triggered LFP: ',tfile]);
        figImage = sprintf('%s%s%s%s',sessions{ii},strcat('STAPlots','\'),cells{jj}(1:end-2),'.fig');
        bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('STAPlots','\'),cells{jj}(1:end-2),'.bmp');
        epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('STAPlots','\'),cells{jj}(1:end-2),'.eps');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');
        saveas(gcf,epsImage,'eps');
    end
end




