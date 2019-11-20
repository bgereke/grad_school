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

function TFRCOUNT(inFile,freqVec,width,powbins)  % width - number of cycles in wavelet (> 5 advisable); order - whitening

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
        cscList = str
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

for ii = 1:numsessions
    disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));
    % Check if subdir for storing images are present. If not, it is
    % created
    dirInfo = dir(sessions{ii});
    found = 0;
    for kk=1:size(dirInfo,1)
        if dirInfo(kk).isdir
            if strcmp(dirInfo(kk).name,strcat('TFRplots','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(sessions{ii},strcat('TFRplots','\'));
    end
    
    % Load data from the .ncs files, make plots, and store them
    for jj=1:numchannels
        disp('Make plots and store them to files');
        disp(sprintf('%s%i',' CSC ',jj, '(indexes number of channels you added)'));
        file = [sessions{ii},channels{jj}];
        [samples,ts,tt, Fs, bv, ir] = loadEEG2(file);
        ch_X = bv*samples;
        [TFR,timeVec,freqVec] = traces2TFR([ch_X ch_X],freqVec,Fs,width);
        [TFRnorm, TFRcounts] = histTFR(TFR,freqVec,powbins);
        TFRcounts = TFRcounts/Fs;
        data.time = timeVec; data.freq = freqVec; data.TFR = TFR;   %save data as struct
        data.TFRnorm = TFRnorm; data.TFRcounts = TFRcounts;
        filename = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\'),channels{jj}(1:end-4),'.mat');
        save(filename,'-struct','data');
        figure(1)
        imagesc(powbins,freqVec,TFRcounts,[0 0.005*max(max(TFRcounts))]); axis xy
        xlabel('normalized power');ylabel('frequency (Hz)');title('Power-Frequency Histogram')
        colorbar
        figImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_powfreqhist_'),channels{jj}(1:end-4),'.fig');
        bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_powfreqhist_'),channels{jj}(1:end-4),'.bmp');
        epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_powfreqhist_'),channels{jj}(1:end-4),'.eps');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');
        saveas(gcf,epsImage,'eps');
        figure(2)
        plot(freqVec,max(TFR,[],2));xlabel('frequency (Hz)');ylabel('max power');
        figImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_maxpow_'),channels{jj}(1:end-4),'.fig');
        bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_maxpow_'),channels{jj}(1:end-4),'.bmp');
        epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_maxpow_'),channels{jj}(1:end-4),'.eps');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');
        saveas(gcf,epsImage,'eps');
        figure(3)
        theta_idx = find(max(TFR(3:7,:),[],2) == max(max(TFR(3:7,:))));
        theta_idx = theta_idx+2;
        v = TFRcounts(theta_idx,1:floor(0.75*(length(powbins)))); % plot out to 75% theta power contour
        contour(powbins,freqVec,TFRcounts,v);
        xlabel('normalized power'); ylabel('frequency (Hz)'); title('power-frequency contour');
        figImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_powfreqcont_'),channels{jj}(1:end-4),'.fig');
        bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_powfreqcont_'),channels{jj}(1:end-4),'.bmp');
        epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_powfreqcont_'),channels{jj}(1:end-4),'.eps');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');
        saveas(gcf,epsImage,'eps');
        [decays] = expfit(TFRcounts,powbins);
        figure(4)
        plot(freqVec,decays);xlabel('frequency (Hz)');ylabel('exponential decay fit');
        figImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_decayfit_'),channels{jj}(1:end-4),'.fig');
        bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_decayfit_'),channels{jj}(1:end-4),'.bmp');
        epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_decayfit_'),channels{jj}(1:end-4),'.eps');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');
        saveas(gcf,epsImage,'eps');
%         figure(5)
%         plot(freqVec,var(TFRnorm,0,2)./mean(TFRnorm,2));xlabel('frequency (Hz)');ylabel('Fano Factor');
%         figImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_fano_'),channels{jj}(1:end-4),'.fig');
%         bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_fano_'),channels{jj}(1:end-4),'.bmp');
%         epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_fano_'),channels{jj}(1:end-4),'.eps');
%         saveas(gcf,figImage,'fig');
%         saveas(gcf,bmpImage,'bmp');
%         saveas(gcf,epsImage,'eps');
        figure(6)
        plot(freqVec,mean(TFRnorm,2)./std(TFRnorm,0,2));xlabel('frequency (Hz)');ylabel('Signal-to-Noise Ratio');
        figImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_snr_'),channels{jj}(1:end-4),'.fig');
        bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_snr_'),channels{jj}(1:end-4),'.bmp');
        epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('TFRplots','\TFR_snr_'),channels{jj}(1:end-4),'.eps');
        saveas(gcf,figImage,'fig');
        saveas(gcf,bmpImage,'bmp');
        saveas(gcf,epsImage,'eps');
    end
end




