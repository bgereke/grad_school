function [sgauto, fgauto, cross] = gamma_xcorr(inFile,lags,sgll,sgl,sgh,sghh,fgll,fgl,fgh,fghh)

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

% Load data from the .ncs files, make plots, and store them
for jj=1%:numchannels
    TFR = [];TP = [];X = [];
    for ii = 1%:numsessions
        xfile = [sessions{ii},channels{jj}];
        [x,~,tt, Fs, bv, ir] = loadEEG2(xfile);
%         x = bv*x;

        SG = abs(hilbert(fftbandpass(x,Fs,sgll,sgl,sgh,sghh)));
        FG = abs(hilbert(fftbandpass(x,Fs,fgll,fgl,fgh,fghh)));
        
        sgauto = acf(SG',SG',lags);
        fgauto = acf(FG',FG',lags);
        crossl = acf(SG',FG',lags);
        crossr = acf(FG',SG',lags);
        cross = [flipud(crossl(2:end));crossr];
    end    
end

