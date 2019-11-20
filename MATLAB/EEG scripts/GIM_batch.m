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

function [GIM] = GIM_batch(inFile,numperms)

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
    mousenum = str(14:15);
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
mousenum = str2num(mousenum);

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
            if strcmp(dirInfo(kk).name,strcat('GIM','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(sessions{ii},strcat('GIM','\'));
    end
        
    load(strcat(sessions{ii},'CS\CS.mat'));
    
    MC = zeros(numchannels,numchannels,length(freqVec));
    pMC = zeros(numchannels,numchannels,length(freqVec),numperms);
    start = round(linspace(310,size(CS,2),numperms));
    
    for jj = 1:numchannels
        for kk = jj:numchannels 
            MC(jj,kk,:) = sum(cs,2);                   
            for p = 1:numperms
                pMC(jj,kk,:,p) = sum(cs-imag(cs)*1i+randn(size(cs))*1i,2);
            end
        end
    end
    
    GIM = zeros(1,length(freqVec));
    MCWPLI = zeros(1,length(freqVec));
    keyboard
%     pGIM = zeros(p,length(freqVec));
    for f = 1:length(freqVec)
        mc = (MC(:,:,f) + tril(MC(:,:,f)',-1)); 
        mc = mc./mean(mean(abs(mc)));
        GIM(f) = 0.5*trace(inv(real(mc))*imag(mc)*inv(real(mc))*imag(mc)');
    end
    
    close all
    figure(1)
    plot(freqVec,GIM,'-k','LineWidth',2);hold on
    plot(freqVec,MCWPLI,'-r','LineWidth',2)
    legend('GIM','MCWPLI')
    xlabel('frequency (Hz)');ylabel('GIM');
    ylim([min([GIM MCWPLI]) max([GIM MCWPLI])])

    figImage = sprintf('%s%s%s%s',sessions{ii},strcat('GIM','\GIM'),'.fig');
    bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('GIM','\GIM'),'.bmp');
    saveas(gcf,figImage,'fig');
    saveas(gcf,bmpImage,'bmp');
    
    clear data wsumImX
end



