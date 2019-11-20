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

function [rpd] = RPD_batch(inFile,freqVec,numbins)

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
            if strcmp(dirInfo(kk).name,strcat('RPD','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(sessions{ii},strcat('RPD','\'));
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
    
    for jj = 1:numchannels
        
        % Load data from the .ncs files, make plots, and store them
        
        disp('Make plots and store them to files');
        xfile = [sessions{ii},channels{jj}];
        [samples,ts,tt, Fs, bv, ir] = loadEEG2(xfile);
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
        
        for kk = jj+1:numchannels
            yfile = [sessions{ii},channels{kk}];
            [samples,ts,tt, Fs, bv, ir] = loadEEG2(yfile);
            y = bv*samples;
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
            [rpd,phasebins] = RPD(x,y,freqVec,6,Fs,numbins);
            %     data.freq = freqVec; data.phasebins = phasebins; data.rpd = rpd;   %save data as struct
            %     data.rpdhi = rpdhi; data.rpdlo = rpdlo;
            %     filename = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'.mat');
            %     save(filename,'-struct','data');
            imagesc(phasebins,freqVec,rpd); axis xy; colormap hot; colorbar
            hold on;
            %     means = rpd*repmat(phasebins(1:end-1)',1,length(phasebins)-1);
            [means] = circ_mean(repmat(phasebins,length(freqVec),1),rpd,2);
            plot(means,freqVec,'k');
            plot(zeros(1,length(freqVec)),freqVec,'k');
            xlabel('relative phase (rads)');ylabel('frequency (Hz)');
            title(strcat('Relative Phase Distribution: ',channels{jj}(1:end-4),'x',channels{kk}(1:end-4)));
            figImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.fig');
            bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'.bmp');
            %     epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'.eps');
            saveas(gcf,figImage,'fig');
            saveas(gcf,bmpImage,'bmp');
            %     figure;
            %     imagesc(phasebins(1:end-1),freqVec,rpd_full); axis xy; colorbar
            %     hold on;
            % %     means = rpd*repmat(phasebins(1:end-1)',1,length(phasebins)-1);
            %     [means] = circ_mean(repmat(phasebins(1:end-1),length(freqVec),1),rpd_full,2);
            %     plot(means,freqVec,'k');
            %     plot(zeros(1,length(freqVec)),freqVec,'k');
            %     xlabel('relative phase (rads)');ylabel('frequency (Hz)');
            %     title(strcat('Relative Phase Distribution: ',channels{1}(1:end-4),'x',channels{2}(1:end-4)));
            %     figImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'full.fig');
            %     bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'full.bmp');
            % %     epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'.eps');
            %     saveas(gcf,figImage,'fig');
            %     saveas(gcf,bmpImage,'bmp');
            %     mpd = max(rpd,[],2);
            %     figure;plot(freqVec,mpd);xlabel('frequency (Hz)');ylabel('rpd max');title(strcat('Relative Phase Distribution Max: ',channels{1}(1:end-4),'x',channels{2}(1:end-4)));
            %     figImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'maxima.fig');
            %     bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'maxima.bmp');
            % %     epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'.eps');
            %     saveas(gcf,figImage,'fig');
            %     saveas(gcf,bmpImage,'bmp');
            [R] = circ_r(repmat(phasebins,length(freqVec),1),rpd,[],2);
            hold off;plot(freqVec,R,'-k');xlabel('frequency (Hz)');ylabel('relative phase vector length');title(strcat('Relative Phase Vector Lengths: ',channels{jj}(1:end-4),'x',channels{kk}(1:end-4)));
            figImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'r.fig');
            bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{jj}(1:end-4),'_',channels{kk}(1:end-4),'r.bmp');
            %     epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'.eps');
            saveas(gcf,figImage,'fig');
            F = getframe(gcf);
            imwrite(F.cdata,bmpImage,'bmp');
            %     saveas(gcf,epsImage,'eps');
            %     close all;figure(2);
            %     imagesc(phasebins(1:end-1),freqVec,rpdhi); axis xy; colorbar
            %     hold on;
            %     means = rpdhi*repmat(phasebins(1:end-1)',1,length(phasebins)-1);
            %     plot(means,freqVec,'k')
            %     plot(zeros(1,length(freqVec)),freqVec,'k');
            %     xlabel('relative phase (rads)');ylabel('frequency (Hz)');
            %     title(strcat('Relative Phase Distribution (z>=2): ',channels{1}(1:end-4),'x',channels{2}(1:end-4)));
            %     figImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'zhi.fig');
            %     bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'zhi.bmp');
            %     epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'zhi.eps');
            %     saveas(gcf,figImage,'fig');
            %     saveas(gcf,bmpImage,'bmp');
            %     saveas(gcf,epsImage,'eps');
            %     close all;figure(3);
            %     imagesc(phasebins(1:end-1),freqVec,rpdlo); axis xy; colorbar
            %     hold on;
            %     means = rpdlo*repmat(phasebins(1:end-1)',1,length(phasebins)-1);
            %     plot(means,freqVec,'k')
            %     plot(zeros(1,length(freqVec)),freqVec,'k');
            %     xlabel('relative phase (rads)');ylabel('frequency (Hz)');
            %     title(strcat('Relative Phase Distribution (z<=0): ',channels{1}(1:end-4),'x',channels{2}(1:end-4)));
            %     figImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'zlo.fig');
            %     bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'zlo.bmp');
            %     epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'zlo.eps');
            %     saveas(gcf,figImage,'fig');
            %     saveas(gcf,bmpImage,'bmp');
            %     saveas(gcf,epsImage,'eps');
            %     close all;figure(4);
            %     imagesc(phasebins(1:end-1),freqVec,rrpdlo); axis xy; colorbar
            %     hold on;
            %     means = rrpdlo*repmat(phasebins(1:end-1)',1,length(phasebins)-1);
            %     plot(means,freqVec,'k')
            %     plot(zeros(1,length(freqVec)),freqVec,'k');
            %     xlabel('relative phase (rads)');ylabel('frequency (Hz)');
            %     title(strcat('Relative Phase Distribution (rz<=0): ',channels{1}(1:end-4),'x',channels{2}(1:end-4)));
            %     figImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'rzlo.fig');
            %     bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'rzlo.bmp');
            %     epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'rzlo.eps');
            %     saveas(gcf,figImage,'fig');
            %     saveas(gcf,bmpImage,'bmp');
            %     saveas(gcf,epsImage,'eps');
            %     close all;figure(5);
            %     imagesc(phasebins(1:end-1),freqVec,rrpdhi); axis xy; colorbar
            %     hold on;
            %     means = rrpdhi*repmat(phasebins(1:end-1)',1,length(phasebins)-1);
            %     plot(means,freqVec,'k')
            %     plot(zeros(1,length(freqVec)),freqVec,'k');
            %     xlabel('relative phase (rads)');ylabel('frequency (Hz)');
            %     title(strcat('Relative Phase Distribution (rz>=2): ',channels{1}(1:end-4),'x',channels{2}(1:end-4)));
            %     figImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'rzhi.fig');
            %     bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'rzhi.bmp');
            %     epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'rzhi.eps');
            %     saveas(gcf,figImage,'fig');
            %     saveas(gcf,bmpImage,'bmp');
            %     saveas(gcf,epsImage,'eps');
            %     close all;figure(6);
            %     imagesc(phasebins(1:end-1),freqVec,irpdlo); axis xy; colorbar
            %     hold on;
            %     means = irpdlo*repmat(phasebins(1:end-1)',1,length(phasebins)-1);
            %     plot(means,freqVec,'k')
            %     plot(zeros(1,length(freqVec)),freqVec,'k');
            %     xlabel('relative phase (rads)');ylabel('frequency (Hz)');
            %     title(strcat('Relative Phase Distribution (iz<=0): ',channels{1}(1:end-4),'x',channels{2}(1:end-4)));
            %     figImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'izlo.fig');
            %     bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'izlo.bmp');
            %     epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'izlo.eps');
            %     saveas(gcf,figImage,'fig');
            %     saveas(gcf,bmpImage,'bmp');
            %     saveas(gcf,epsImage,'eps');
            %     close all;figure(7);
            %     imagesc(phasebins(1:end-1),freqVec,irpdhi); axis xy; colorbar
            %     hold on;
            %     means = irpdhi*repmat(phasebins(1:end-1)',1,length(phasebins)-1);
            %     plot(means,freqVec,'k')
            %     plot(zeros(1,length(freqVec)),freqVec,'k');
            %     xlabel('relative phase (rads)');ylabel('frequency (Hz)');
            %     title(strcat('Relative Phase Distribution (iz>=2): ',channels{1}(1:end-4),'x',channels{2}(1:end-4)));
            %     figImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'izhi.fig');
            %     bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'izhi.bmp');
            %     epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'izhi.eps');
            %     saveas(gcf,figImage,'fig');
            %     saveas(gcf,bmpImage,'bmp');
            %     saveas(gcf,epsImage,'eps');
            %     close all;figure(8);
            %     imagesc(phasebins(1:end-1),freqVec,pfr); axis xy; colorbar
            %     xlabel('relative phase (rads)');ylabel('frequency (Hz)');
            %     title(strcat('z-scored cross spectral apmlitude: ',channels{1}(1:end-4),'x',channels{2}(1:end-4)));
            %     figImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'pfr.fig');
            %     bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'pfr.bmp');
            %     epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'pfr.eps');
            %     saveas(gcf,figImage,'fig');
            %     saveas(gcf,bmpImage,'bmp');
            %     saveas(gcf,epsImage,'eps');
            %     close all;figure(9);
            %     imagesc(phasebins(1:end-1),freqVec,rpfr); axis xy; colorbar
            %     xlabel('relative phase (rads)');ylabel('frequency (Hz)');
            %     title(strcat('z-scored cross spectral apmlitude (real): ',channels{1}(1:end-4),'x',channels{2}(1:end-4)));
            %     figImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'rpfr.fig');
            %     bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'rpfr.bmp');
            %     epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'rpfr.eps');
            %     saveas(gcf,figImage,'fig');
            %     saveas(gcf,bmpImage,'bmp');
            %     saveas(gcf,epsImage,'eps');
            %     close all;figure(10);
            %     imagesc(phasebins(1:end-1),freqVec,ipfr); axis xy; colorbar
            %     xlabel('relative phase (rads)');ylabel('frequency (Hz)');
            %     title(strcat('z-scored cross spectral apmlitude (imaginary): ',channels{1}(1:end-4),'x',channels{2}(1:end-4)));
            %     figImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'ipfr.fig');
            %     bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'ipfr.bmp');
            %     epsImage = sprintf('%s%s%s%s',sessions{ii},strcat('RPD','\'),channels{1}(1:end-4),'_',channels{2}(1:end-4),'ipfr.eps');
            %     saveas(gcf,figImage,'fig');
            %     saveas(gcf,bmpImage,'bmp');
            %     saveas(gcf,epsImage,'eps');
        end
    end
end




