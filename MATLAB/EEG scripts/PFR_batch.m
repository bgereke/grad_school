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

function [cycle] = PFR_batch(inFile,freqVec,phasebins)

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
        str(3:7) = [];
        str(1) = 'G';
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
        str(3:7) = [];
        str(1) = 'G';
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

% Set the field selection for reading the video files. 1 = Add parameter, 0 = skip
% parameter
fieldSelection(1) = 1; % Timestamps
fieldSelection(2) = 1; % Extracted X
fieldSelection(3) = 1; % Extracted Y
fieldSelection(4) = 0; % Extracted Angel
fieldSelection(5) = 1; % Targets
fieldSelection(6) = 0; % Points
% Do we return header 1 = Yes, 0 = No.
extractHeader = 0;
% 5 different extraction modes, see help file for Nlx2MatVt
extractMode = 1; % Extract all data

V = [];
D = 200;

for ii = 1%:numsessions
    disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));
    % Check if subdir for storing images are present. If not, it is
    % created
    dirInfo = dir(sessions{ii});
    found = 0;
    for kk=1:size(dirInfo,1)
        if dirInfo(kk).isdir
            if strcmp(dirInfo(kk).name,strcat('PFR','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(sessions{ii},strcat('PFR','\'));
    end
   
    
    % Load data from the .ncs files, make plots, and store them
    for jj=1:numchannels
        disp('Make plots and store them to files');
        disp(sprintf('%s%i',' CSC ',jj,' of ',numchannels));
        file = [sessions{ii},channels{jj}];
        [samples,ts,tt, Fs, bv, ir] = loadEEG2(file);
        ch_X = bv*samples;  
        clear samples

        [pfr_mean,cycle] = PFR(ch_X,freqVec,6,Fs,phasebins);
        
        figure(1);
        col = max(max(abs(pfr_mean)));
        surf(phasebins,freqVec,zeros(size(pfr_mean)),pfr_mean);view([0 90]);axis tight;axis square;shading interp;
        caxis([-col col]);colorbar
        set(gca,'Ytick',0:10:250,'yticklabel',0:10:250,'YScale','linear','XScale','linear');
        xlabel('theta phase');ylabel('Frequency (Hz)');
        title(strcat('Theta Phase Frequency Representation: ',channels{jj}));hold on
        plot(phasebins,100*(cycle-min(cycle))/max(cycle-min(cycle))+10,'-k','LineWidth',2);
        text(0, 120, num2str(max(cycle)),'FontSize',12);hold off
        bmpImage = sprintf('%s%s%s%s',strcat(sessions{ii},'PFR','\'),channels{jj}(1:end-4),'.bmp');
        saveas(gcf,bmpImage,'bmp');        
    end
end




