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
% C:\Data\cscList.txt
% C:\Data\Begin 1
% C:\Data\Begin 2
% C:\Data\Begin 3
% C:\Data\Begin 4
% and so on ...
%
% 'CSCList.txt' contains a list of the Neuralynx .csc files to be analyzed.
% All plots will be stored to both bmp and eps imagefiles to a subdirectory in
% the data folder called CFC_plots.

function gammawindows(inFile,s1,s2,f1,f2)

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
            if strcmp(dirInfo(kk).name,strcat('GammaWindowPlots','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(sessions{ii},strcat('GammaWindowPlots','\'));
    end
    
    % Load data from the .ncs files, make plots, and store them
    for jj=1:numchannels
        disp('Make plots and store them to files');
        disp(sprintf('%s%i',' CSC ',jj, ' of ', numchannels));
        file = [sessions{ii},channels{jj}];
        [samples,ts,tt, Fs, bv, ir] = loadEEG2(file);
        ch_X = bv*samples;
        % slow gamma
        figure(1);
        [peak_ind_unique,non_overlap_peak_ind_unique,start2,stop2,non_overlap_gamma_windows_CA1] = code_for_slow_fast_gamma_indices_windows3_CA1only(ch_X, s1, s2, Fs);
        data_s.s1 = s1; data_s.s2 = s2; data_s.swindows = non_overlap_gamma_windows_CA1;   %save data as struct
        filename = sprintf('%s%s%s%s',sessions{ii},strcat('GammaWindowPlots','\slow_'),channels{jj}(1:end-4),'.mat');
        save(filename,'-struct','data_s');
        xlabel('time (ms)');ylabel('uV');
        title(['Slow Gamma Window ',num2str(s1),'-',num2str(s2),' Hz: ',file]);
        figImage = sprintf('%s%s%s%s',sessions{ii},strcat('GammaWindowPlots','\slow_'),channels{jj}(1:end-4),'.fig');
        bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('GammaWindowPlots','\slow_'),channels{jj}(1:end-4),'.bmp');
        epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('GammaWindowPlots','\slow_'),channels{jj}(1:end-4),'.eps');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');
        saveas(gcf,epsImage,'eps');
        % fast gamma
        figure(1);
        [peak_ind_unique,non_overlap_peak_ind_unique,start2,stop2,non_overlap_gamma_windows_CA1] = code_for_slow_fast_gamma_indices_windows3_CA1only(ch_X, f1, f2, Fs);
        data_f.f1 = f1; data_f.f2 = s2; data_f.fwindows = non_overlap_gamma_windows_CA1;   %save data as struct
        filename = sprintf('%s%s%s%s',sessions{ii},strcat('GammaWindowPlots','\fast_'),channels{jj}(1:end-4),'.mat');
        save(filename,'-struct','data_f');
        xlabel('time (ms)');ylabel('uV');
        title(['Fast Gamma Window ',num2str(f1),'-',num2str(f2),' Hz: ',file]);
        figImage = sprintf('%s%s%s%s',sessions{ii},strcat('GammaWindowPlots','\fast_'),channels{jj}(1:end-4),'.fig');
        bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('GammaWindowPlots','\fast_'),channels{jj}(1:end-4),'.bmp');
        epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('GammaWindowPlots','\fast_'),channels{jj}(1:end-4),'.eps');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');
        saveas(gcf,epsImage,'eps');
    end
end




