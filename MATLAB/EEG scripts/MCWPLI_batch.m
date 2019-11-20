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

function [MCWPLI] = MCWPLI_batch(inFile,freqVec,width,minutes,numperms)

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
            if strcmp(dirInfo(kk).name,strcat('MCWPLI','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(sessions{ii},strcat('MCWPLI','\'));
    end
        
    Wx = cell(numchannels,1);
    for jj = 1:numchannels
        % Load data from the .ncs files, make plots, and store them
        xfile = [sessions{ii},channels{jj}];
        [samples,~,tt, Fs, bv, ~] = loadEEG2(xfile);
        x = bv*samples;
        %set recording channel against ground and then against average reference
%         if mousenum ~= 47
%             xcscnum = channels{jj}(4:end-4); % get  tt# corresponding to avgsjj
%             xcscnum = str2num(xcscnum);
%             xttnum = ceil(xcscnum/4);
%             if ~strcmp(refs{xttnum},'G')
%                 xrfile = [sessions{ii},refs{xttnum}];
%                 [samples,~,tt, Fs, bv, ~] = loadEEG2(xrfile);
%                 x = x + bv*samples;
%             end
%         else
%             xcscnum = channels{jj}(4:end-4); % get  tt# corresponding to avgsjj
%             xcscnum = str2num(xcscnum);
%             xttnum = ceil(xcscnum/4);
%             if xttnum == 6 || xttnum == 10
%                 xrfile = [sessions{ii},refs{xttnum}];
%                 [samples,~,~,~, bv,~] = loadEEG2(xrfile);
%                 x = x - bv*samples;
%             end
%         end
        ss = 50;
        x = x(tt-min(tt)<=minutes*60);
        tt = tt(tt-min(tt)<=minutes*60);
        tt = tt(1:ss:end);
        Wx{jj} = traces2Wx(x,freqVec,Fs,width);
        Wx{jj} = Wx{jj}(:,1:ss:end);
    end
    
    MC = zeros(numchannels,numchannels,length(freqVec));
    iMC = zeros(numchannels,numchannels,length(freqVec));
    aiMC = zeros(numchannels,numchannels,length(freqVec));
%     pMC = zeros(numchannels,numchannels,length(freqVec),numperms);
%     awsumImX = zeros(size(Wx{1}));
%     wsumImX = zeros(size(Wx{1}));
%     pwsumImX = zeros(size(Wx{1}));
    
    for jj = 1:numchannels
        if jj > 1
           Wx{jj-1} = []; %clear used wavelet transforms
        end
        for kk = jj:numchannels
            
            cs = Wx{jj}.*conj(Wx{kk});

%             if jj ~= kk
%                 wpli = sum(imag(cs),2)./sum(abs(imag(cs)),2);
%                 %             dwpli = sqrt(sum((imag(Wx{jj}.*conj(Wx{kk}))).^2,2));
%                 %             WPLI = [WPLI abs(wpli)];
%                 wsumImX = wsumImX+imag(cs).*repmat(wpli,1,size(Wx{jj},2));  
% %                 add = abs(imag(cs).*repmat(wpli,1,size(Wx{jj},2))).*sign(rand(size(cs))-0.5);
% %                 pwsumImX = pwsumImX + add.*repmat(sign(mean(add,2)),1,size(cs,2));
%             end
            
            MC(jj,kk,:) = sum(cs,2); 
            iMC(jj,kk,:) = abs(sum(imag(cs),2));
            aiMC(jj,kk,:) = sum(abs(imag(cs)),2);            
%             
%             for p = 1:numperms
%                 pMC(jj,kk,:,p) = sum(cs-imag(cs)*1i+randn(size(cs))*1i,2);
%             end
            
        end
    end
    
%     wMCWPLI = abs(sum(wsumImX,2))./sum(abs(wsumImX),2);
%     pwMCWPLI = abs(sum(pwsumImX,2))./sum(abs(pwsumImX),2);
    
    GIM = zeros(1,length(freqVec));
    MCWPLI = zeros(1,length(freqVec));
%     pGIM = zeros(p,length(freqVec));
    for f = 1:length(freqVec)
        mc = (MC(:,:,f) + tril(MC(:,:,f)',-1)); 
        mc = mc./mean(mean(abs(mc)));
        GIM(f) = 0.5*trace(inv(real(mc))*imag(mc)*inv(real(mc))*imag(mc)');
        
        num = (iMC(:,:,f) + tril(iMC(:,:,f)',-1));
        den = (aiMC(:,:,f) + tril(aiMC(:,:,f)',-1));
        e = eig(inv(den)*num)
        MCWPLI(f) = sum(abs(e))/numchannels;
%         e = eig(inv(den)*num);
%         MCWPLI(f) = max(max(e));
%         e = eig(inv(real(mc))*imag(mc)*inv(real(mc))*imag(mc)');
%         GIM(f) = sqrt(0.5*sum(e.^2));
%         for p = 1:numperms
%             pmc = (pMC(:,:,f,p) + tril(pMC(:,:,f,p)',-1));
%             pmc = pmc./mean(mean(abs(mc)));
%             pGIM(p,f) = 0.5*trace(inv(real(pmc))*imag(pmc)*inv(real(pmc))*imag(pmc)');
%         end
    end

    %     wMCWPLId = abs(sum(wsumImXd,2))./sum(abs(wsumImXd),2);
    %     awMCWPLI = abs(sum(wsumImX,2))./sum(awsumImX,2);
    %     sumWPLI = sum(WPLI,2)/max(sum(WPLI,2));
    %     wWPLI = sum(WPLI.^2,2)./sum(WPLI,2);
    %     maxWPLI = max(WPLI,[],2);
    
    close all
    figure(1)
%     subplot(1,2,1)
%     [ci] = prctile(pGIM,[2.5 97.5]);
%     shadedplot(freqVec, ci(2,:), ci(1,:),[1 0 0],[1 0.25 0.25]); alpha(0.5);hold on
%     plot(freqVec,mean(pGIM),'-r')
    plot(freqVec,GIM,'-k','LineWidth',2);hold on
    plot(freqVec,MCWPLI,'-r','LineWidth',2)
    legend('GIM','MCWPLI')
    xlabel('frequency (Hz)');ylabel('GIM');
    ylim([min([GIM MCWPLI]) max([GIM MCWPLI])])
%     
%     subplot(1,2,2)
%     plot(freqVec,(GIM-mean(pGIM))./std(pGIM),'-k','LineWidth',2)
%     xlabel('frequency (Hz)');ylabel('GIM (z-score from null)');
%     
%     figure(2)
%     subplot(1,2,1)
%     plot(freqVec,pwMCWPLI,'-r','LineWidth',2)
%     plot(freqVec,wMCWPLI,'-k','LineWidth',2)
%     xlabel('frequency (Hz)');ylabel('wMCWPLI');

    %     plot(freqVec,wMCWPLId,'-r','LineWidth',2);
    %     plot(freqVec,awMCWPLI,'-r','LineWidth',2);
    %     plot(freqVec,sumWPLI,'-b','LineWidth',2);
    %     plot(freqVec,wWPLI,'-g','LineWidth',2);
    %     plot(freqVec,maxWPLI,'-m','LineWidth',2);
    %     legend('MCWPLI','aMCWPLI','sumWPLI','wWPLI','maxWPLI');
%     data.freq = freqVec; data.wmcwpli = wMCWPLI; data.MCImX = wsumImX; data.tt = tt; %data.awmcwpli = awMCWPLI;   %save data as struct
%     filename = sprintf('%s%s%s%s',sessions{ii},strcat('MCWPLI','\MCWPLI'),'.mat');
%     save(filename,'-struct','data');
    figImage = sprintf('%s%s%s%s',sessions{ii},strcat('MCWPLI','\MCWPLI'),'.fig');
    bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('MCWPLI','\MCWPLI'),'.bmp');
    saveas(gcf,figImage,'fig');
    saveas(gcf,bmpImage,'bmp');
    
    clear data wsumImX
end



