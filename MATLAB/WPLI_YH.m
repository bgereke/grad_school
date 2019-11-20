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
% C:\Data\CSCList.txt (first CSC is DG, the second is CA3)
% C:\Data\Begin 1
% C:\Data\Begin 2
% C:\Data\Begin 3
% C:\Data\Begin 4
% and so on ...
%
% 'CSCList.txt' contains a list of the Neuralynx .csc files to be analyzed.
% All plots will be stored to both bmp and eps imagefiles to a subdirectory in
% the data folder called CFC_plots.
% (-pi:10/360*2*pi:pi)=phasebins

function [WPLI] = WPLI_YH(inFile,freqVec,width)

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

bootWPLI = cell(numsessions,1);
WPLI = cell(numsessions,1);

for ii = 1:numsessions
    disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));
    % Check if subdir for storing images are present. If not, it is
    % created
    dirInfo = dir(sessions{ii});
    found = 0;
    for kk=1:size(dirInfo,1)
        if dirInfo(kk).isdir
            if strcmp(dirInfo(kk).name,strcat('WPLIplots','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(sessions{ii},strcat('WPLIplots','\'));
    end
    
    % Load data from the .ncs files, make plots, and store them
    disp('Make plots and store them to files'); 
    xfile = [sessions{ii},channels{1}];
    [samples,ts,tt, Fs, bv, ir] = loadEEG2(xfile);
    x = bv*samples;
    yfile = [sessions{ii},channels{2}];
    [samples,ts,tt, Fs, bv, ir] = loadEEG2(yfile);
    y = bv*samples;
%         freqVec = 1:2:100;
%         phasebins = -pi:10/360*2*pi:pi;

    clear samples

    [ImX] = traces2ImX(x,y,freqVec,Fs,width);
    ImX = ImX(:,1:19:end);
    outsum   = nansum(ImX,2);      % compute the sum; 
    outsumW  = nansum(abs(ImX),2); % normalization of the WPLI
    outssq   = nansum(ImX.^2,2);
    WPLI{ii} = abs(outsum.^2 - outssq)./(outsumW.^2 - outssq); % do the pairwise thing in a handy way
%     WPLI = WPLI.^(0.5);

    %block bootstrap
    tic
    Bstar = opt_block_length_REV_dec07(ImX'); %optimal block lengths for each frequency
    toc
    %circular
    tic
    cImX  = [ImX ImX(:,1:max(max(ceil(Bstar))))]; %wrap ImX around on itself
    numboots = 10000;
    bootWPLI{ii} = zeros(length(freqVec),numboots);
    for b = 1:numboots
        for f = 1:length(freqVec)
            %create blocks
            bstart = ceil(size(ImX,2)*rand(1,ceil(size(ImX,2)/ceil(Bstar(2,f)))));
            blocks = repmat(bstart,ceil(Bstar(2,f)),1) + repmat([0:ceil(Bstar(2,f))-1]',1,size(bstart,2));
            bootImX = cImX(f,blocks(:)); %sample blocks
            bootImX = bootImX(1:size(ImX,2)); %remove extra samples
            outsum   = nansum(bootImX,2);      % compute the sum; 
            outsumW  = nansum(abs(bootImX),2); % normalization of the WPLI
            outssq   = nansum(bootImX.^2,2);
            bootWPLI{ii}(f,b) = abs(outsum.^2 - outssq)./(outsumW.^2 - outssq); %debiased squared WPLI
        end
        if ~mod(b,100)
            disp(sprintf('%s%s',num2str(round(b/numboots*100)),'%'));
        end
    end
    toc    
    
end

WPLI = (WPLI{1}+WPLI{2}+WPLI{3})/3;

upper = prctile((bootWPLI{1}+bootWPLI{2}+bootWPLI{3})/3,97.5,2);
lower = prctile((bootWPLI{1}+bootWPLI{2}+bootWPLI{3})/3,2.5,2);

plot(freqVec,WPLI,'k');hold on;plot(freqVec,upper,'--k');plot(freqVec,lower,'--k')
title(strcat('WPLI: ',channels{1}(1:end-4),'x',channels{2}(1:end-4)));

data.freqvec = freqVec; data.wpli = WPLI; data.upper = upper; data.lower = lower;  %save data as struct
filename = sprintf('%s%s%s%s',sessions{ii},strcat('WPLIplots','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'.mat');
save(filename,'-struct','data');
close all
figure(1)
plot(freqVec,WPLI,'-k'); xlabel('Frequency'); ylabel('Debiased Squared WPLI');hold on
plot(freqVec,upper,'--k');plot(freqVec,lower,'--k');
title(strcat('WPLI: ',channels{1}(1:end-4),'x',channels{2}(1:end-4)));
figImage = sprintf('%s%s%s%s',sessions{ii},strcat('WPLIplots','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'.fig');
bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('WPLIplots','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'.bmp');
epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('WPLIplots','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'.eps');
saveas(gcf,figImage,'fig');
saveas(gcf,bmpImage,'bmp');
saveas(gcf,epsImage,'epsc');
clear WPLI data

