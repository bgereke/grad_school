%convert CSC's from .ncs to .dat for use in spikedetekt 

function [] = ncs2dat(inFile)

img_text = 'on';

fid = fopen(inFile,'r');
if fid == -1
    msgbox('Could not open the input file! Make sure the filname and path are correct.','ERROR');
end

parent = cd;

% Get sessions and csc-file list from input file
fid = fopen(inFile,'r');
ii = -1;     
numsessions = 0;
while ~feof(fid)
    str = fgetl(fid);
    if ii == -1
        cscList = str;
    else
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
       if ~strcmp(str(end),'s')
            str = strcat(str,'.ncs');
        end
       channels(jj) = {str};
       jj = jj+1;
end
numchannels = jj-1;
cscid = fclose('all');

%save .dat file for each shank
numshanks = str2num(channels{numchannels}(2));
% Fs = 32000;
% [b,a] = butter(3, [500 0.95*0.5*Fs]/(Fs/2), 'bandpass'); %3rd order highpass butterworth filter
for s = 1:numshanks
    shank = [];
    disp(sprintf('%s%s','Reading data for shank: ',num2str(s)));
    %contatenate sessions
    for ii = 1:numsessions
        xfile = [sessions{ii},channels{1}];
        [~,~,tt, ~, ~, ~] = loadEEG2(xfile);
        cscs = zeros(numchannels,length(tt)+1); clear tt
        for jj = 1:numchannels
            xfile = [sessions{ii},channels{jj}];
            [samples,~,~, ~, bv, ~] = loadEEG2(xfile);
            cscs(jj,2:end) = bv*samples;%filter(b,a,bv*samples); 
            cscs(jj,1) = str2num(channels{jj}(2));clear samples
        end
        %median reference to remove movement artifacts
        shidx = find(cscs(:,1)==s);
        catshank = zeros(length(shidx),size(cscs,2)-1);
        med = median(cscs(:,2:end));
        for jj = 1:length(shidx)
            catshank(jj,:) = cscs(shidx(jj),2:end)-med;
        end
        clear cscs
        shank = [shank catshank]; clear catshank med
    end

    filename = strcat('S',num2str(s),'CSCs.dat');
    fid = fopen(filename,'w');
    fwrite(fid,shank,'int16');
    fclose(fid);
end