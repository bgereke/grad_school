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

function [PFR] = PFR_ImX_batch(inFile,freqVec,phasebins)

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
            if strcmp(dirInfo(kk).name,strcat('PFR_ImX','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(sessions{ii},strcat('PFR_ImX','\'));
    end
    
     % Make Average Reference
%     avref = 0;
%     for jj=1:numrefs
%         file = [sessions{ii},avgs{jj}];
%         [samples,ts,tt, Fs, bv, ir] = loadEEG2(file);
%         ch_X = bv*samples;
%         cscnum = avgs{jj}(4:end-4); % get  tt# corresponding to avgsjj 
%         cscnum = str2num(cscnum);   
%         ttnum = ceil(cscnum/4);
%         if ~strcmp(refs{ttnum},'G')
%             rfile = [sessions{ii},refs{ttnum}];
%             [samples,ts,tt, Fs, bv, ir] = loadEEG2(rfile);
%             ch_X = ch_X + bv*samples;
%         end
%         avref = avref + ch_X/numrefs;
%         clear samples ch_X
%     end
    
    % Load data from the .ncs files, make plots, and store them
    
    disp('Make plots and store them to files'); 
    xfile = [sessions{ii},channels{1}];
    [samples,ts,tt, Fs, bv, ir] = loadEEG2(xfile);
    x = bv*samples;
    yfile = [sessions{ii},channels{2}];
    [samples,ts,tt, Fs, bv, ir] = loadEEG2(yfile);
    y = bv*samples;
    %set recording channel against ground and then against average reference
    xcscnum = channels{1}(4:end-4); % get  tt# corresponding to avgsjj
    xcscnum = str2num(xcscnum);
    xttnum = ceil(xcscnum/4);
    ycscnum = channels{1}(4:end-4); % get  tt# corresponding to avgsjj
    ycscnum = str2num(ycscnum);
    yttnum = ceil(ycscnum/4);
    if ~strcmp(refs{xttnum},'G')
        xrfile = [sessions{ii},refs{xttnum}];
        [samples,ts,tt, Fs, bv, ir] = loadEEG2(xrfile);
        x = x + bv*samples;
    end
    if ~strcmp(refs{yttnum},'G')
        yrfile = [sessions{ii},refs{yttnum}];
        [samples,ts,tt, Fs, bv, ir] = loadEEG2(yrfile);
        y = y + bv*samples;
    end
%     x = x - avref;
%     y = y - avref;
    clear samples
    [PFR] = PFR_ImX(x,y,freqVec,6,Fs,phasebins);
%     PFR = PFR.^(0.5);
    data.freq = freqVec; data.phasebins = phasebins; data.pfr = PFR;   %save data as struct
    filename = sprintf('%s%s%s%s',sessions{ii},strcat('PFR_ImX','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'.mat');
    save(filename,'-struct','data');
    figure;
    imagesc(phasebins(1:end-1),freqVec,PFR); axis xy; colorbar
    xlabel('phase (rads)');ylabel('frequency (Hz)');
    title(strcat('Average Z-scored Theta Phase Triggered ImX: ',channels{1}(1:end-4),'x',channels{2}(1:end-4)));
    figImage = sprintf('%s%s%s%s',sessions{ii},strcat('PFR_ImX','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'.fig');
    bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('PFR_ImX','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'.bmp');
    epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('PFR_ImX','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'.eps');
    saveas(gcf,figImage,'fig');
    saveas(gcf,bmpImage,'bmp');
    saveas(gcf,epsImage,'eps');
    
end




