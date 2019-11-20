% WPLI('inFile.txt')
%
% WPLI generates time frequency representaion of the weighted phase lag index
% plots it and stores it to images and .mat files.
% Multiple sessions will be read from the CSC's specififed in
% 'CSCList.txt'. 'TTList.txt' is necessary for spike data scripts that use
% the same 'infile.txt'.
%
% The input file must be on the following format.
%
% C:\Data\TTList.txt
% C:\Data\CSCList.txt
% C:\Data\refList.txt
% C:\Data\ch4avg.txt
% C:\Data\Begin 1
% C:\Data\Begin 2
% C:\Data\Begin 3
% C:\Data\Begin 4
% and so on ...
%
% 'CSCList.txt' contains a list of the Neuralynx .csc files to be analyzed.
% All plots will be stored to both bmp and eps imagefiles to a subdirectory in
% the data folder called CFC_plots.

function crossWPLI(inFile,freqVec,width,cycles,minutes)  % width - number of cycles in wavelet (> 5 advisable); order - whitening


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
    for jj = 1:numchannels
        disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));
        % Check if subdir for storing images are present. If not, it is
        % created
        dirInfo = dir(sessions{ii});
        found = 0;
        for kk=1:size(dirInfo,1)
            if dirInfo(kk).isdir
                if strcmp(dirInfo(kk).name,strcat('crossWPLI','\'))
                    found = 1;
                end
            end
        end
        if found==0
            mkdir(sessions{ii},strcat('crossWPLI','\'));
        end
        
        %     % Make Average Reference
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
        
        filex = [sessions{ii},channels{jj}];
        [samples,ts,tt, Fs, bv, ir] = loadEEG2(filex);
        x = bv*samples;
        %set recording channel against ground and then against average reference
        xcscnum = channels{jj}(4:end-4); % get  tt# corresponding to avgsjj
        xcscnum = str2num(xcscnum);
        xttnum = ceil(xcscnum/4);
        if ~strcmp(refs{xttnum},'G')
            xrfile = [sessions{ii},refs{xttnum}];
            [samples,ts,tt, Fs, bv, ir] = loadEEG2(xrfile);
            x = x + bv*samples;
        end
        x = x(tt-min(tt)<=minutes*60);
        Wx = traces2Wx(x,freqVec,Fs,width);
        
        L = cycles/min(freqVec);
        dl = 1/(20*max(freqVec));
        tlags = -L:dl:L;
        clags = -cycles:0.001:cycles;
        xWPLI_t  = zeros(length(freqVec),length(tlags));
        xWPLI_p  = zeros(length(freqVec),length(clags));
        
        for kk = jj+1:numchannels
            filey = [sessions{ii},channels{kk}];
            [samples,ts,tt, Fs, bv, ir] = loadEEG2(filey);
            y = bv*samples;
            ycscnum = channels{kk}(4:end-4); % get  tt# corresponding to avgsjj
            ycscnum = str2num(ycscnum);
            yttnum = ceil(ycscnum/4);
            if ~strcmp(refs{yttnum},'G')
                yrfile = [sessions{ii},refs{yttnum}];
                [samples,ts,tt, Fs, bv, ir] = loadEEG2(yrfile);
                y = y + bv*samples;
            end
            y = y(tt-min(tt)<=minutes*60);
            Wy = traces2Wx(y,freqVec,Fs,width);
            Wy = conj(Wy);
            %     x = x - avref;
            %     y = y - avref;
            clear samples   
            
            for l = 1:length(tlags)
                
                shiftidx = round(tlags(l)*Fs);    
                if shiftidx < 0
                    ImX = imag(Wx(:,1:end+shiftidx).*Wy(:,-shiftidx+1:end));
                elseif shiftidx > 0
                    ImX = imag(Wx(:,shiftidx+1:end).*Wy(:,1:end-shiftidx));
                else
                    ImX = imag(Wx.*Wy);
                end
                outsum   = nansum(ImX,2);      % compute the sum;
                outsumW  = nansum(abs(ImX),2); % normalization of the WPLI
                xWPLI_t(:,l) = (outsum./outsumW).^2;
            end
            for f = 1:length(freqVec)
                for c = 1:length(clags)
                    if min(abs(tlags-clags(c)/freqVec(f)))==0
                        xWPLI_p(f,c) = xWPLI_t(f,abs(tlags-clags(c)/freqVec(f))==0);
                    else
                        idx1 = find(tlags-clags(c)/freqVec(f)>0,1,'first');
                        p1 = abs(tlags(idx1)-clags(c)/freqVec(f));
                        idx2 = find(tlags-clags(c)/freqVec(f)<0,1,'last');
                        p2 = abs(tlags(idx2)-clags(c)/freqVec(f));
                        xWPLI_p(f,c) = (xWPLI_t(f,idx1)/p1+xWPLI_t(f,idx2)/p2)/(1/p1+1/p2);
                    end
                end                
            end
            
            imagesc(clags,freqVec,xWPLI_p);axis xy;colormap hot;colorbar
            xlabel('lag (cycles)');ylabel('frequency (Hz)');
            title(strcat('cross-WPLI: ',channels{jj}(1:end-4),'x',channels{kk}(1:end-4)));
            data.freqvec = freqVec; data.xwpli = xWPLI_p;  %save data as struct
            filename = sprintf('%s%s%s%s',sessions{ii},strcat('crossWPLI','\'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.mat');
            save(filename,'-struct','data');
            figImage = sprintf('%s%s%s%s',sessions{ii},strcat('crossWPLI','\'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.fig');
            bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('crossWPLI','\'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.bmp');
            saveas(gcf,figImage,'fig');
            saveas(gcf,bmpImage,'bmp');
            clear xWPLI_t xWPLI_p data
        end
    end
end




