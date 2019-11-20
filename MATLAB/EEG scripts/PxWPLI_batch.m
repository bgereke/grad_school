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

function [PxWPLI] = PxWPLI_batch(inFile,freqVec,numbins)

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
    foundm = 0;
    foundz = 0;
    founds = 0;
    for kk=1:size(dirInfo,1)
        if dirInfo(kk).isdir
            if strcmp(dirInfo(kk).name,strcat('PxWPLI','\'))
                found = 1;
            end
            if strcmp(dirInfo(kk).name,strcat('mPxWPLI','\'))
                foundm = 1;
            end
            if strcmp(dirInfo(kk).name,strcat('zPxWPLI','\'))
                foundz = 1;
            end
            if strcmp(dirInfo(kk).name,strcat('sPxWPLI','\'))
                founds = 1;
            end
        end
    end
    if found==0
        mkdir(sessions{ii},strcat('PxWPLI','\'));
    end
    if foundm==0
        mkdir(sessions{ii},strcat('mPxWPLI','\'));
    end
    if foundz==0
        mkdir(sessions{ii},strcat('zPxWPLI','\'));
    end
    if founds==0
        mkdir(sessions{ii},strcat('sPxWPLI','\'));
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
    
    for jj = 1:numchannels-1
        % Load data from the .ncs files, make plots, and store them
        
        disp('Make plots and store them to files');
        xfile = [sessions{ii},channels{jj}];
        [samples,ts,tt, Fs, bv, ir] = loadEEG2(xfile);
        x = bv*samples;
        %set recording channel against ground and then against average reference
        xcscnum = channels{1}(4:end-4); % get  tt# corresponding to avgsjj
        xcscnum = str2num(xcscnum);
        xttnum = ceil(xcscnum/4);
        if ~strcmp(refs{xttnum},'G')
            xrfile = [sessions{ii},refs{xttnum}];
            [samples,ts,tt, Fs, bv, ir] = loadEEG2(xrfile);
            x = x + bv*samples;
        end
        
        for kk = jj+1:numchannels
            yfile = [sessions{ii},channels{kk}];
            [samples,ts,tt, Fs, bv, ir] = loadEEG2(yfile);
            y = bv*samples;
            %set recording channel against ground and then against average referenc
            ycscnum = channels{kk}(4:end-4); % get  tt# corresponding to avgsjj
            ycscnum = str2num(ycscnum);
            yttnum = ceil(ycscnum/4);            
            if ~strcmp(refs{yttnum},'G')
                yrfile = [sessions{ii},refs{yttnum}];
                [samples,ts,tt, Fs, bv, ir] = loadEEG2(yrfile);
                y = y + bv*samples;
            end
            %     x = x - avref;
            %     y = y - avref;
            clear samples
            abins = linspace(-2,2,numbins);
            [PxWPLI] = getPxWPLI(x,y,freqVec,5,Fs,abins);            
            data.freq = freqVec; data.abins = abins; data.PxWPLI = PxWPLI;   %save data as struct
            filename = sprintf('%s%s%s%s',sessions{ii},strcat('PxWPLI','\'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.mat');
            save(filename,'-struct','data');
            figure(1);%subplot(1,2,1)
            imagesc(freqVec,abins,PxWPLI(:,:,1)); axis xy square; colormap hot;colorbar
            xlabel('frequency (Hz)');ylabel('z-scored power');
            title(strcat('Average WPLI vs Power: ',channels{jj}(1:end-4),channels{kk}(1:end-4)));
%             subplot(1,2,2)
%             imagesc(freqVec,abins,PxWPLI(:,:,2)); axis xy square; colormap hot;colorbar
%             xlabel('frequency (Hz)');ylabel('z-scored power');
%             title(strcat('Average WPLI vs Power: ',channels{kk}(1:end-4)));
            figImage = sprintf('%s%s%s%s',sessions{ii},strcat('PxWPLI','\'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.fig');
            bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('PxWPLI','\'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.bmp');
            saveas(gcf,figImage,'fig');
            saveas(gcf,bmpImage,'bmp');
%             data.freq = freqVec; data.abins = abins; data.mpfr = mPFR; data.pfr = [];   %save data as struct
%             filename = sprintf('%s%s%s%s',sessions{ii},strcat('mPxWPLI','\'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.mat');
%             save(filename,'-struct','data');
%             figure(2);%subplot(1,2,1)
%             xcm = max(max(mPxWPLI(:,freqVec>=50,1)));
% %             ycm = max(max(mPxWPLI(:,freqVec>=50,2)));
%             imagesc(freqVec,abins,mPxWPLI(:,:,1),[-xcm xcm]); axis xy square; colormap jet;colorbar
%             xlabel('frequency (Hz)');ylabel('z-scored power');
%             title(strcat('Change in WPLI vs Power: ',channels{jj}(1:end-4)));
% %             subplot(1,2,2)
% %             imagesc(freqVec,abins,mPxWPLI(:,:,2),[-ycm ycm]); axis xy square; colormap jet;colorbar
% %             xlabel('frequency (Hz)');ylabel('z-scored power');
% %             title(strcat('Change in WPLI vs Power: ',channels{kk}(1:end-4)));
%             figImage = sprintf('%s%s%s%s',sessions{ii},strcat('mPxWPLI','\'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.fig');
%             bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('mPxWPLI','\'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.bmp');
%             saveas(gcf,figImage,'fig');
%             saveas(gcf,bmpImage,'bmp');
%             data.freq = freqVec; data.abins = abins; data.zpfr = zPFR; data.mpfr = [];   %save data as struct
%             filename = sprintf('%s%s%s%s',sessions{ii},strcat('zPxWPLI','\'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.mat');
%             save(filename,'-struct','data');
%             figure(3);subplot(1,2,1)
%             xcz = prctile(reshape(zPxWPLI(:,:,1),size(zPxWPLI(:,:,1),1)*size(zPxWPLI(:,:,1),2),1),95);
%             ycz = prctile(reshape(zPxWPLI(:,:,2),size(zPxWPLI(:,:,2),1)*size(zPxWPLI(:,:,2),2),1),95);
%             imagesc(freqVec,abins,zPxWPLI(:,:,1),[-xcz,xcz]); axis xy square; colormap jet;colorbar
%             xlabel('frequency (Hz)');ylabel('z-scored power');
%             title(strcat('Z-Scored Change in WPLI vs Power: ',channels{jj}(1:end-4)));
%             subplot(1,2,2)
%             imagesc(freqVec,abins,zPxWPLI(:,:,2),[-ycz,ycz]); axis xy square; colormap jet;colorbar
%             xlabel('frequency (Hz)');ylabel('z-scored power');
%             title(strcat('Z-Scored Change in WPLI vs Power: ',channels{kk}(1:end-4)));
%             figImage = sprintf('%s%s%s%s',sessions{ii},strcat('zPxWPLI','\'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.fig');
%             bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('zPxWPLI','\'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.bmp');
%             saveas(gcf,figImage,'fig');
%             saveas(gcf,bmpImage,'bmp');
%             data.freq = freqVec; data.abins = abins; data.spfr = sPFR; data.zpfr = [];   %save data as struct
%             filename = sprintf('%s%s%s%s',sessions{ii},strcat('sPxWPLI','\'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.mat');
%             save(filename,'-struct','data');
%             figure(4);subplot(1,2,1)
%             imagesc(freqVec,abins,sPxWPLI(:,:,1),[-xcm xcm]); axis xy square; colormap jet;colorbar
%             xlabel('frequency (Hz)');ylabel('z-scored power');
%             title(strcat('Significant Change in WPLI vs Power: ',channels{jj}(1:end-4)));
%             subplot(1,2,2)
%             imagesc(freqVec,abins,sPxWPLI(:,:,2),[-ycm ycm]); axis xy square; colormap jet;colorbar
%             xlabel('frequency (Hz)');ylabel('z-scored power');
%             title(strcat('Significant Change in WPLI vs Power: ',channels{kk}(1:end-4)));
%             figImage = sprintf('%s%s%s%s',sessions{ii},strcat('sPxWPLI','\'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.fig');
%             bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('sPxWPLI','\'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.bmp');
%             saveas(gcf,figImage,'fig');
%             saveas(gcf,bmpImage,'bmp');
        end
    end

end




