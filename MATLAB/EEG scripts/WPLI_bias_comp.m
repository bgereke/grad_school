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

function [MCWPLI] = WPLI_bias_comp(inFile,freqVec,width,minutes)

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
%     found = 0;
%     for kk=1:size(dirInfo,1)
%         if dirInfo(kk).isdir
%             if strcmp(dirInfo(kk).name,strcat('MCWPLI','\'))
%                 found = 1;
%             end
%         end
%     end
%     if found==0
%         mkdir(sessions{ii},strcat('MCWPLI','\'));
%     end
    
   
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
            x = x(tt-min(tt)<=minutes*60);
            Wx{jj} = traces2Wx(x,freqVec,Fs,width); 
            wsumX = zeros(size(Wx{jj}));
            csumX = zeros(size(Wx{jj}));
            sumwpli = zeros(size(Wx{jj},1),1);
            sumcv = zeros(size(Wx{jj},1),1);
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
                x = x(tt-min(tt)<=minutes*60);
                Wx{kk} = traces2Wx(x,freqVec,Fs,width);
                clear samples
            end 
            
            wpli = sum(imag(Wx{jj}.*conj(Wx{kk})),2)./sum(abs(imag(Wx{jj}.*conj(Wx{kk}))),2);
            cv = mean(imag(Wx{jj}.*conj(Wx{kk})),2)./std(abs(imag(Wx{jj}.*conj(Wx{kk}))),0,2);
            sumwpli = sumwpli+abs(wpli);
            sumcv = sumcv+abs(cv);
            wsumX = wsumX+Wx{jj}.*conj(Wx{kk}).*repmat(wpli,1,size(Wx{jj},2));
            csumX = csumX+Wx{jj}.*conj(Wx{kk}).*repmat(cv,1,size(Wx{jj},2));
        end
    end
    wsumX = wsumX./repmat(sumwpli,1,size(wsumX,2));
    csumX = csumX./repmat(sumcv,1,size(csumX,2));  
    
    reps = 1000;
    smax = 2000;
    ds = 2:10:smax;
    basicwpli = zeros(reps,length(ds));
    debiasedwpli = zeros(reps,length(ds));
    mywpli = zeros(reps,length(ds));
    mypli = zeros(reps,length(ds));
    mycv = zeros(reps,length(ds));
    
    for s = 1:length(ds)      
        for r = 1:reps
            start = randsample(size(wsumX,2)-ds(s),1);
            sidx = start:start+ds(s)-1;
            basicwpli(r,s) = abs(sum(imag(wsumX(20,sidx))))/sum(abs(imag(wsumX(20,sidx))));
            debiasedwpli(r,s) = (sum(imag(wsumX(20,sidx))).^2-sum(imag(wsumX(20,sidx)).^2))/(sum(abs(imag(wsumX(20,sidx)))).^2-sum(imag(wsumX(20,sidx)).^2));
            mywpli(r,s) = sum(imag(wsumX(20,sidx)))/sum(abs(imag(wsumX(20,sidx))));
            mypli(r,s) = sum(sign(imag(wsumX(20,sidx))))/sum(sign(abs(imag(wsumX(20,sidx)))));
            mycv(r,s) = mean(imag(csumX(20,sidx)))./std(imag(csumX(20,sidx)));
        end        
    end    
%     mycv = mycv/std(abs(imag(csumX(20,sidx))),0,2);
    
    keyboard
    plot(ds,mean(basicwpli),'k');hold on
    plot(ds,sqrt(mean(debiasedwpli)),'b');
    plot(ds,mean(mywpli),'r');
    plot(ds,mean(mypli),'m');
    plot(ds,mean(mycv),'g');
    
    
%     figure(1);hold off
%     plot(freqVec,wMCWPLI,'-k','LineWidth',2);xlabel('frequency (Hz)');ylabel('MCWPLI');
%     hold on
%     plot(freqVec,awMCWPLI,'-r','LineWidth',2);
%     plot(freqVec,sumWPLI,'-b','LineWidth',2);  
%     plot(freqVec,wWPLI,'-g','LineWidth',2);
%     plot(freqVec,maxWPLI,'-m','LineWidth',2);
%     legend('MCWPLI','aMCWPLI','sumWPLI','wWPLI','maxWPLI');
%     data.freq = freqVec; data.wmcwpli = wMCWPLI;data.awmcwpli = awMCWPLI;   %save data as struct
%     filename = sprintf('%s%s%s%s',sessions{ii},strcat('MCWPLI','\MCWPLI'),'.mat');
%     save(filename,'-struct','data');
%     figImage = sprintf('%s%s%s%s',sessions{ii},strcat('MCWPLI','\MCWPLI'),'.fig');
%     bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('MCWPLI','\MCWPLI'),'.bmp');
%     saveas(gcf,figImage,'fig');
%     saveas(gcf,bmpImage,'bmp');
end




