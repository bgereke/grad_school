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

function [bestCSCs] = bestCSCs_batch(inFile,freqVec,width,minutes)

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
            if strcmp(dirInfo(kk).name,strcat('bestCSCs','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(sessions{ii},strcat('bestCSCs','\'));
    end
    
    numpairs = sum(1:numchannels-1);
    drange = zeros(length(freqVec),numpairs);
    chnums = zeros(2,numpairs);
    bestCSCs = zeros(length(freqVec),2);
    pairnum = 1;
    Wx = cell(numchannels,1);
    Phi = zeros(length(freqVec),1);
    Plo = zeros(length(freqVec),1);
    
    for jj = 1:numchannels-1
        % Load data from the .ncs files, make plots, and store them
        
        disp('Make plots and store them to files');
        if jj == 1
            xfile = [sessions{ii},channels{jj}];
            [samples,~,tt, Fs, bv, ~] = loadEEG2(xfile);
            x = bv*samples;
            %set recording channel against ground and then against average reference
            xcscnum = channels{1}(4:end-4); % get  tt# corresponding to avgsjj
            xcscnum = str2num(xcscnum);
            xttnum = ceil(xcscnum/4);
            if ~strcmp(refs{xttnum},'G')
                xrfile = [sessions{ii},refs{xttnum}];
                [samples,~,tt, Fs, bv, ~] = loadEEG2(xrfile);
                x = x + bv*samples;
            end
            x = x(tt-min(tt)<=minutes*60);
            Wx{jj} = traces2Wx(x,freqVec,Fs,width);            
        end
        
        for kk = jj+1:numchannels
            if jj == 1
                xfile = [sessions{ii},channels{kk}];
                [samples,~,tt, Fs, bv, ~] = loadEEG2(xfile);
                x = bv*samples;
                %set recording channel against ground and then against average referenc
                xcscnum = channels{kk}(4:end-4); % get  tt# corresponding to avgsjj
                xcscnum = str2num(xcscnum);
                xttnum = ceil(xcscnum/4);
                if ~strcmp(refs{xttnum},'G')
                    xrfile = [sessions{ii},refs{xttnum}];
                    [samples,~,tt, Fs, bv, ~] = loadEEG2(xrfile);
                    x = x + bv*samples;
                end
                x = x(tt-min(tt)<=minutes*60);
                Wx{kk} = traces2Wx(x,freqVec,Fs,width);
                clear samples
            end            
            ImX = imag(Wx{jj}.*conj(Wx{kk}));
%             tfr = (2*abs(Wx{kk})/Fs).^2;
%             med = median(tfr,2);
%             for f = 1:length(freqVec)
%                 outsum = sum(ImX(f,tfr(f,:)>med(f)),2);
%                 outsumW = sum(abs(ImX(f,tfr(f,:)>med(f))),2);
%                 Phi(f) = (outsum./outsumW).^2;
%                 outsum = sum(ImX(f,tfr(f,:)<med(f)),2);
%                 outsumW = sum(abs(ImX(f,tfr(f,:)<med(f))),2);
%                 Plo(f) = (outsum./outsumW).^2;
%             end
%             drange(:,pairnum) = Phi-Plo;
            drange(:,pairnum) = (sum(ImX,2)./sum(abs(ImX),2)).^2;
            chnums(:,pairnum) = [jj;kk];
            pairnum = pairnum+1;    
        end
    end
    
    [~,bidx] = max(drange,[],2);
    bestnums = chnums(:,bidx)';    
    for f = 1:length(freqVec)
        for c = 1:2
            cscnum = channels{bestnums(f,c)}(4:end-4);
            bestCSCs(f,c) = str2num(cscnum);
        end
    end
    
    figure(1)
    for c = 1:size(drange,2)
        plot(freqVec,drange(:,c)','-k'); xlabel('frequency (Hz)'); ylabel('max of joint WPLI')
        hold on
    end
    md = max(drange,[],2);
    plot(freqVec,md,'-r','LineWidth',2)
    hold off
    data.freq = freqVec; data.drange = drange; data.bestCSCs = bestCSCs; data.chnums = chnums;  %save data as struct
    filename = sprintf('%s%s%s%s',sessions{ii},strcat('bestCSCs','\best'),'.mat');
    save(filename,'-struct','data');
    figImage = sprintf('%s%s%s%s',sessions{ii},strcat('bestCSCs','\maxjointWPLI'),'.fig');
    bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('bestCSCs','\maxjointWPLI'),'.bmp');
    saveas(gcf,figImage,'fig');
    saveas(gcf,bmpImage,'bmp');
end




