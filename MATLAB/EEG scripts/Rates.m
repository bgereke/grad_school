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

function [Rates] = Rates(inFile)

img_text = 'on';

fid = fopen(inFile,'r');
if fid == -1
    msgbox('Could not open the input file! Make sure the filname and path are correct.','ERROR');
end

% Get sessions and csc-file list from input file
fid = fopen(inFile,'r');
ii = 0;     
while ~feof(fid)
    str = fgetl(fid);
    if ii == 0
        ttList = str;
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

Rates = zeros(numsessions,numcells);

for ii = 1:numsessions
    disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));
    
    % Load data from the .tt files, make plots, and store them
    for jj=1:numcells
        disp('Make plots and store them to files');
        disp(sprintf('%s%i',' cell ',jj));
            
        dirFiles = dir(sessions{ii});
       
        if isempty(strmatch(char(cells{jj}),char(dirFiles.name),'exact'))
            continue
        end
         
        tfile = [sessions{ii},cells{jj}];
        [TS] = loadSpikes(tfile);              
        
        % Get position data
        % Set the field selection for reading the video files. 1 = Add parameter, 0 = skip
        % parameter
        fieldSelection(1) = 1; % Timestamps
        fieldSelection(2) = 1; % Extracted X
        fieldSelection(3) = 1; % Extracted Y
        fieldSelection(4) = 0; % Extracted Angel
        fieldSelection(5) = 0; % Targets
        fieldSelection(6) = 0; % Points
        % Do we return header 1 = Yes, 0 = No.
        extractHeader = 0;
        % 5 different extraction modes, see help file for Nlx2MatVt
        extractMode = 1; % Extract all data
        file = strcat(sessions{ii},'vt1.nvt');
        [t, x, y] = Nlx2MatVT(file,fieldSelection,extractHeader,extractMode);
        % Convert timestamps to seconds
        t = t/1000000;
        t = max(t)-min(t); %length of session
        
        Rates(ii,jj) = length(TS)/t;
    end
end





