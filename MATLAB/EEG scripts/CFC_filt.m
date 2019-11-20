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

function [CFC_filt] = CFC_filt(inFile,freqVec)

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
    if ii == 0
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

% read the file names from the csc-file list
cscid = fopen(cscList,'r');
jj = 1;
while ~feof(cscid)
       str = fgetl(cscid);
       channels(jj) = {str};
       jj = jj+1;
end
numchannels = jj-1;

for ii = 1:numsessions
    disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));
    % Check if subdir for storing images are present. If not, it is
    % created
    dirInfo = dir(sessions{ii});
    found = 0;
    for kk=1:size(dirInfo,1)
        if dirInfo(kk).isdir
            if strcmp(dirInfo(kk).name,strcat('CFC_filt_plots','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(sessions{ii},strcat('CFC_filt_plots','\'));
    end
    
    % Load data from the .ncs files, make plots, and store them
    for jj=1:numchannels
        disp('Make plots and store them to files');
        disp(sprintf('%s%i',' CSC ',jj, ' of ',numchannels));
        file = [sessions{ii},channels{jj}];
        [samples,ts,tt, Fs, bv, ir] = loadEEG2(file);
        ch_X = bv*samples;
        
        [pass_times, pass_idx] = thetaDeltaFilter(ch_X, tt, Fs,[], 2.5);
        
        if size(pass_idx) == [0,0]
            continue
        end
        
        [CFC_filt,freqVec1,freqVec2] = crossfreqCoh(ch_X(pass_idx(1,1):pass_idx(1,2)),ch_X(pass_idx(1,1):pass_idx(1,2)),freqVec,Fs);
        
        for ww = 2:size(pass_idx,1)
            [CFC,freqVec1,freqVec2] = crossfreqCoh(ch_X(pass_idx(ww,1):pass_idx(ww,2)),ch_X(pass_idx(ww,1):pass_idx(ww,2)),freqVec,Fs);
            CFC_filt = CFC_filt + (pass_idx(ww,2)-pass_idx(ww,1))*CFC;
        end
        
        CFC_filt = CFC_filt/sum(pass_idx(ww,2)-pass_idx(ww,1));
        
        data.freq1 = freqVec1; data.freq2 = freqVec2; data.crf = CFC_filt;   %save data as struct
        filename = sprintf('%s%s%s%s',sessions{ii},strcat('CFC_filt_plots','\'),channels{jj}(1:end-4),'.mat');
        save(filename,'-struct','data');
        figure(1);
        imagesc(freqVec1,freqVec2,CFC_filt); axis xy; xlim([0 max(freqVec)]);
        xlabel('Frequency_p_h_a_s_e (Hz)');ylabel('Frequency_a_m_p_l_i_t_u_d_e (Hz)');
        title(strcat('Cross Frequency Coherence: ',file));
        figImage = sprintf('%s%s%s%s',sessions{ii},strcat('CFC_filt_plots','\'),channels{jj}(1:end-4),'.fig');
        bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('CFC_filt_plots','\'),channels{jj}(1:end-4),'.bmp');
        epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('CFC_filt_plots','\'),channels{jj}(1:end-4),'.eps');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');
        saveas(gcf,epsImage,'eps');
    end
end




