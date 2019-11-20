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

function [MCWPLI] = traces2MCWPLI(inFile,freqVec,width)

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
            if strcmp(dirInfo(kk).name,strcat('MCWPLI_series','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(sessions{ii},strcat('MCWPLI_series','\'));
    end
    
   
    Wx = cell(numchannels,1);
    
    for jj = 1:numchannels-1
        % Load data from the .ncs files, make plots, and store them
        
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
            Wx{jj} = traces2Wx(x,freqVec,Fs,width); 
            wsumImX = zeros(size(Wx{jj}));
        else
            Wx{jj-1} = [];
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
                Wx{kk} = traces2Wx(x,freqVec,Fs,width);
                clear samples
            end 
            
            wpli = sum(imag(Wx{jj}.*conj(Wx{kk})),2)./sum(abs(imag(Wx{jj}.*conj(Wx{kk}))),2);
            wsumImX = wsumImX+imag(Wx{jj}.*conj(Wx{kk})).*repmat(wpli,1,size(Wx{jj},2));
        end
    end
    
    MCWPLI = zeros(size(wsumImX));
    for j=1:length(freqVec)
        MCWPLI(j,:) = wavepli(freqVec(j),wsumImX(j,:),Fs,width);
    end
    data.freq = freqVec; data.wmcwpli = MCWPLI;   %save data as struct
    filename = sprintf('%s%s%s%s',sessions{ii},strcat('MCWPLI_series','\MCWPLI'),'.mat');
    save(filename,'-struct','data');
    clear MCWPLI data wsumImX
end

function wpli = wavepli(f,ImX,Fs,width)
% function wpli = wpli(f,Wx,Wy,Fs,width)
%
% Return a vector containing the wpli as a
% function of time for frequency f. The wpli
% is calculated by integrating over a gaussian window around time t
% set by width. 
% ImX: imaginary part of the cross spectrum (frequency x time)
% Fs: sampling frequency
% width : width of integtration window (>= 5 suggested).
%
% Brian Gereke, February 2015

dt = 1/Fs;
sf = f/width;
st = 1/(2*pi*sf);

t=-3.5*st:dt:3.5*st;
g = gaussian(f,t,width);
os = convfft(ImX,g);
osW = convfft(abs(ImX),g);
% ossq = convfft(ImX.^2,g);
os = os(ceil(length(g)/2):length(os)-floor(length(g)/2));
osW = osW(ceil(length(g)/2):length(osW)-floor(length(g)/2));
% ossq = ossq(ceil(length(g)/2):length(ossq)-floor(length(g)/2));
wpli = os./(osW); 

function g = gaussian(f,t,width)
% function y = morlet(f,t,width)
% 
% Morlet's wavelet for frequency f and time t. 
% The wavelet will be normalized so the total energy is 1.
% width defines the ``width'' of the wavelet. 
% A value >= 5 is suggested.
%
% Ref: Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997)
%
% See also: PHASEGRAM, PHASEVEC, WAVEGRAM, ENERGY 
%
% Ole Jensen, August 1998 

sf = f/width;
st = 1/(2*pi*sf);
A = (st*sqrt(pi))^(-0.5);

g = A*exp(-t.^2/(2*st^2));




