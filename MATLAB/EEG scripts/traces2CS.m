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

function [CS] = traces2CS(inFile,freqVec,width)

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
        ttList = str;
    elseif ii == 0
%         str(1) = 'E';
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
%         str(1) = 'E';
        sessions(numsessions) = {str};
    end
    ii = ii+1;
end

% read the file names from the tt-file list
% ttid = fopen(ttList,'r');
% jj = 1;
% 
% while ~feof(ttid)
%     str = fgetl(ttid);
%     cells(jj) = {str};
%     jj = jj+1;
% end
% numcells = jj-1;

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

% pdim = numchannels;
% % read the file names of the references for the channels to be used for
% % averaging
% refid = fopen(refList,'r');
% jj = 1;
% while ~feof(refid)
%     str = fgetl(refid);
%     refs(jj) = {str};
%     jj = jj+1;
% end
% 
% % read the file names of the channels to be used for averaging
% avgid = fopen(ch4avg,'r');
% jj = 1;
% while ~feof(avgid)
%     str = fgetl(avgid);
%     avgs(jj) = {str};
%     jj = jj+1;
% end
% numrefs = jj-1;

% Set the field selection for reading the video files. 1 = Add parameter, 0 = skip
% parameter
fieldSelection(1) = 1; % Timestamps
fieldSelection(2) = 0; % Extracted X
fieldSelection(3) = 0; % Extracted Y
fieldSelection(4) = 0; % Extracted Angel
fieldSelection(5) = 1; % Targets
fieldSelection(6) = 0; % Points
% Do we return header 1 = Yes, 0 = No.
extractHeader = 0;
% 5 different extraction modes, see help file for Nlx2MatVt
extractMode = 1; % Extract all data

D = 96.5; %track diameter in cm
t = cell(numsessions,1);
tt = cell(numsessions,1);
phase = cell(numsessions,1);
vel = cell(numsessions,1);
acc = cell(numsessions,1);

for ii = 1:numsessions
    disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));
    
    % Check if subdir for storing images are present. If not, it is
    % created
    dirInfo = dir(sessions{ii});
    found = 0;
    for kk=1:size(dirInfo,1)
        if dirInfo(kk).isdir
            if strcmp(dirInfo(kk).name,strcat('CS','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(sessions{ii},strcat('CS','\'));
    end
    
    %compute running speed
    % Get position data
    file = strcat(sessions{ii},'vt1.nvt');
    [t{ii}, targets] = Nlx2MatVT(file,fieldSelection,extractHeader,extractMode);
    % Convert timestamps to seconds
    t{ii} = t{ii}'/1000000;
    % Decode the target data
    [dTargets,tracking] = decodeTargets(targets);
    % Exctract position data from the target data
    [frontX,frontY,backX,backY] = extractPosition(dTargets,tracking);
    % Smooth the exctracted position samples using a moving mean filter
    % Set missing samples to NaN. Necessary when using the mean filter.
    fidx = frontX==0&frontY==0; bidx = backX==0&backY==0;
    frontX(fidx) = NaN;frontY(fidx) = NaN;
    backX(bidx) = NaN;backY(bidx) = NaN;
    %some days might only have one recorded LED (via Dylan effect)
    if sum(isnan(frontX))/length(frontX) > 0.5 || sum(isnan(backX))/length(backX)>0.5
        % Set the field selection for reading the video files. 1 = Add parameter, 0 = skip
        % parameter
        fieldSelect(1) = 1; % Timestamps
        fieldSelect(2) = 1; % Extracted X
        fieldSelect(3) = 1; % Extracted Y
        fieldSelect(4) = 0; % Extracted Angel
        fieldSelect(5) = 0; % Targets
        fieldSelect(6) = 0; % Points
        % Do we return header 1 = Yes, 0 = No.
        extractHead = 0;
        % 5 different extraction modes, see help file for Nlx2MatVt
        extractMod = 1; % Extract all data
        [t{ii}, frontX, frontY] = Nlx2MatVT(file,fieldSelect,extractHead,extractMod);
        t{ii} = t{ii}'/1000000;
        ind = find(frontX == 0);
        frontX(ind) = NaN;
        frontY(ind) = NaN;
        backX = frontX;
        backY = frontY;
    end
    
    %scale and center trajectory
    nb = 200;
    xymap = zeros(nb,nb);
    v = 0.05;
    mx = linspace(min(backX),max(backX),nb);
    my = linspace(min(backY),max(backY),nb);
    
    for xx = 1:length(mx)
        dx = mx(xx)-backX;
        for yy = 1:length(my)
            dy = my(yy)-backY;
            xymap(xx,yy) = nansum(v^2/sqrt(2*pi)*exp(-0.5*(dx.*dx*v^2+dy.*dy*v^2)));
        end
    end
    regmax = imregionalmax(xymap,4);
    [xx,yy] = find(regmax);
    [z, a, b, alpha] = fitellipse([mx(xx);my(yy)]); %'linear' to speed up
    %Translate
    frontX = frontX-z(1);backX = backX-z(1);
    frontY = frontY-z(2);backY = backY-z(2);
    %Rotate
    Q = [cos(-alpha), -sin(-alpha); sin(-alpha) cos(-alpha)];
    ftemp = Q*[frontX;frontY];
    btemp = Q*[backX;backY];
    %Scale
    ftemp(1,:) = 0.5*D/a*ftemp(1,:);
    ftemp(2,:) = 0.5*D/b*ftemp(2,:);
    btemp(1,:) = 0.5*D/a*btemp(1,:);
    btemp(2,:) = 0.5*D/b*btemp(2,:);
    %Rotate back to orginal orientation
    Q = [cos(alpha), -sin(alpha); sin(alpha) cos(alpha)];
    ftemp = Q*ftemp;
    btemp = Q*btemp;
    frontX = ftemp(1,:);frontY = ftemp(2,:);
    backX = btemp(1,:);backY = btemp(2,:);
    %Fill missing position samples using polar interpolation
    [fphase,fr] = interporPos(frontX,frontY,2,29.97);
    [bphase,br] = interporPos(backX,backY,2,29.97);
    indfb = isnan(fphase) & isnan(bphase);
    t{ii}(indfb) = [];fphase(indfb) = [];bphase(indfb) = [];
    phase{ii} = circ_mean([fphase;bphase])';
    %smooth and get derivatives
    %     phase{ii} = medfilt1(unwrap(phase{ii}),5);
%     [phase{ii},~] = localfit_ad(t{ii},unwrap(phase{ii}),2,0.04,10);
%     %     phase{ii} = medfilt1(unwrap(phase{ii}),5);
%     [~,vel{ii}] = localfit_ad(t{ii},phase{ii},2,0.175,5);
%     [~,acc{ii}] = localfit_ad(t{ii},abs(vel{ii}),2,0.175,5);
%     t{ii}(isnan(acc{ii})) = [];phase{ii}(isnan(acc{ii})) = [];
%     vel{ii}(isnan(acc{ii})) = [];acc{ii}(isnan(acc{ii})) = [];
%     %     vel{ii} = medfilt1(vel{ii},5);
%     phase{ii} = wraptopi(phase{ii});
%     phase{ii} = reshape(phase{ii},length(phase{ii}),1);
%     vel{ii} = reshape(vel{ii},length(vel{ii}),1);
%     acc{ii} = reshape(acc{ii},length(acc{ii}),1);
%     vel{ii} = vel{ii}*D/2;
    mint = min(t{ii}); maxt = max(t{ii});
    th = 0.035;
    
%     %load spike times/positions/running speeds for all cells
%     TS = [];
%     for jj=1:numcells
%         tfile = [sessions{ii},cells{jj}];
%         [ts] = loadSpikes(tfile);
%         ts = ts(ts>=mint & ts<=maxt);
%         [sptimes] = spikeTimes(ts,t{ii});
%         ts = ts(min(abs(sptimes(:,[2 4])),[],2)<=th);
%         sptimes = sptimes(min(abs(sptimes(:,[2 4])),[],2)<=th,:);
%         TS = [TS;jj*ones(length(ts),1) ts];
%         clear ts spikePos sptimes w
%     end
    
    %load wavelet transforms and get spike phases
    xfile = [sessions{ii},channels{1}];
    [~,~,tt{ii,1}, ~, ~, ~] = loadEEG2(xfile);
%     Wx = zeros(length(freqVec),length(t{ii}),numchannels);
%      Wx = zeros(length(freqVec),length(tt{ii}),numchannels);
    for jj = 1%:numchannels
        % Load data from the .ncs files, make plots, and store them
        xfile = [sessions{ii},channels{jj}];
        [samples,~,~, Fs, bv, ~] = loadEEG2(xfile);
        x = -bv*samples; clear samples
        [theta_phase,wave_phase,sym] = thetaphase(x,Fs);
%         sym(:,1) = tt{ii}(sym(:,1));
%         theta_stuff = [tt{ii}', slowBP',fastBP'];
%     if ii==2,keyboard,end
        kappa = 50;
        pbins = linspace(-pi,pi,60);
        delta = angle(exp(1i*repmat(theta_phase',1,length(pbins))).*conj(exp(1i*repmat(pbins,length(theta_phase),1))));
        W = exp(kappa*cos(delta));
        den = ones(length(freqVec),length(theta_phase))*exp(kappa*cos(delta));
        cycle = x'*W./den;
        cycle = cycle(1,:);

%         if pbins(cycle==min(cycle)) > 0
%             x = -x;
%             [theta_phase,wave_phase,~,~] = thetaphase(x,Fs);
%         end

%         x = std(x)*randn(size(x));
        [wx] = traces2Wx(x,freqVec,Fs,width,'Morlet','area');
  

        if jj == 1
            ttidx = zeros(length(t{ii}),1);
            for it = 1:length(t{ii})
                [~,idx] = min((tt{ii}-t{ii}(it)).^2);
                ttidx(it) = idx;
            end
            
%             TS = TS(TS(:,2)>=min(tt{ii}) & TS(:,2)<=max(tt{ii}),:);
% %             PS = zeros(size(TS,1),length(freqVec));
% %             PS = zeros(size(TS,1),numchannels*length(freqVec));
%             sptimes = spikeTimes(TS(:,2),tt{ii}');
%             w = [sptimes(:,4)./(sptimes(:,4)-sptimes(:,2)) abs(sptimes(:,2))./(sptimes(:,4)-sptimes(:,2))];
        end

        Wx = abs(wx(:,ttidx));
        tt{ii} = tt{ii}(ttidx);
%         Wx = abs(wx);

        phase{ii} = wraptopi(pchip(t{ii},unwrap(phase{ii}),tt{ii}));
        
%         for f = 1:length(freqVec)
%             %             PS(:,(jj-1)*length(freqVec)+f) = exp(1i*circ_mean([angle(wx(f,sptimes(:,1)));angle(wx(f,sptimes(:,3)))],w'));
%             PS(:,(jj-1)*length(freqVec)+f) = circ_mean([angle(wx(f,sptimes(:,1)));angle(wx(f,sptimes(:,3)))],w');
%         end
%         PS = [circ_mean([theta_phase(sptimes(:,1));theta_phase(sptimes(:,3))],w');...
%             circ_mean([wave_phase(sptimes(:,1));wave_phase(sptimes(:,3))],w')];
        theta_phase = theta_phase(ttidx);
        wave_phase = wave_phase(ttidx);
%         theta_phase = theta_phase;
%         wave_phase = wave_phase;
    end
    
%     plot(freqVec,mean(Wx,2));xlabel('frequency (Hz)');ylabel('amplitude (uV)');hold on
%     pause(0.5)
%     clear wx 
%     clear sptimes w
    
    %     %do pca and cross-spectal stuff by frequency to save on memory
    %     TS = [TS zeros(size(TS,1),2*length(freqVec))]; %cell id's, spike times, and phases wrt phasePP
    %     PICS = zeros(length(freqVec),length(tt{ii})); %first principle component of the squared imaginary cross-spectra
    %     CICS = zeros(length(freqVec),length(tt{ii})); %first canonical component of the imaginary cross-spectra
    %     phasePP = zeros(length(freqVec),length(tt{ii})); %first principal component of the complex-valued energy
    % %     dPICS = zeros(length(freqVec),length(tt{ii})); %d -> time derivative
    % %     dCICS = zeros(length(freqVec),length(tt{ii}));
    % %     dPP = zeros(length(freqVec),length(tt{ii}));
    %     pccoeff = zeros(length(freqVec),numchannels);
    %     CS = zeros(length(tt{ii}),numchannels,numchannels); %cross-spectra
    %     for f = 1:length(freqVec)
    %
    %         %get phasePP, dPP, and spike phases
    % %         pcw = pca(squeeze(Wx(f,:,:)),'Centered',false);
    %         [pcw,scores] = pca(squeeze(Wx(f,:,:)));
    %         phasePP(f,:) = scores(:,1);
    %         pccoeff(f,:) = pcw(:,1);
    % %         dPP(f,:) = gradient(abs(phasePP(f,:)),1/Fs);
    %         phaseSpike = PS(:,f:length(freqVec):end)*pcw(:,1);
    %         TS(:,2+f) = real(phaseSpike);
    %         TS(:,2+length(freqVec)+f) = imag(phaseSpike);clear phaseSpike
    %
    %         %compute cross spectra
    %         CS = zeros(length(tt{ii}),numchannels,numchannels);
    %         for jj = 1:numchannels
    %             CS(:,:,jj) = bsxfun(@times,Wx(f,:,:),conj(Wx(f,:,jj)));
    %         end
    %
    %         %get CICS and dCICS
    %         MCf = squeeze(sum(CS,1));
    %         N = size(CS,1);
    %         MCf = MCf/sum(N);
    %         mc = MCf - tril(MCf,-1) + tril(MCf',-1);
    %         re = real(mc)^-0.5;
    %         rho = re*imag(mc)*re;
    %         [V,~] = svd(rho*rho');
    %         alpha = V(:,1); beta = V(:,2);
    %         CICS(f,:) = imag(squeeze(sum(bsxfun(@times,permute(CS,[2 3 1]),(alpha'*re)'),1))'*re*beta); %vectorize the shit out of it
    %         CICS(f,:) = CICS(f,:)*sign(mean(CICS(f,:))); %consider ttidx
    % %         dCICS(f,:) = gradient(CICS(f,:),1/Fs);
    %
    %         %get PICS and dPICS (pca on upper triangle)
    %         numidx = length(find(triu(squeeze(CS(1,:,:)),1)));
    %         [I,J] = ind2sub4up(1:numidx);
    %         pc=pca(imag(reshape(CS(:,sub2ind(size(squeeze(CS(1,:,:))),I,J)),size(CS,1),numidx)).^2);
    %         PICS(f,:) = (imag(reshape(CS(:,sub2ind(size(squeeze(CS(1,:,:))),I,J)),size(CS,1),numidx)).^2)*pc(:,1);
    % %         dPICS(f,:) = gradient(PICS(f,:),1/Fs);
    %
    %         clear CS
    %     end
    %     keyboard
    % %     phasePP = phasePP(:,ttidx);%dPP = dPP(:,ttidx);
    % %     CICS = CICS(:,ttidx); %dCICS = dCICS(:,ttidx);
    % %     PICS = PICS(:,ttidx); %dPICS = dPICS(:,ttidx);
    %     close all
    %     subplot(3,1,1)
    %     plot(freqVec,mean(CICS,2),'-k','LineWidth',2);
    %     xlabel('Frequency (Hz)');ylabel('CICS'); axis square
    %     subplot(3,1,2)
    %     plot(freqVec,log(mean(PICS,2)),'-k','LineWidth',2); axis square
    %     xlabel('Frequency (Hz)')
    %     ylabel('log PICS')
    %     subplot(3,1,3)
    %     plot(freqVec,log(mean(abs(phasePP).^2,2)),'-k','LineWidth',2); axis square
    %     xlabel('Frequency (Hz)')
    %     ylabel('log PP')
    %
    %     figImage = sprintf('%s%s%s%s',sessions{ii},strcat('CS','\GIM'),'.fig');
    %     bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('CS','\GIM'),'.bmp');
    %     saveas(gcf,figImage,'fig');
    %     saveas(gcf,bmpImage,'bmp');

%     %downsample and create lagged maxima   
%     wx = abs(wx);
%     WX = zeros(size(wx));
%     for f = 1:length(freqVec)
%        [~,idx] = findpeaks(wx(f,:),'MinPeakProminence',0.01*mean(wx(:)));
% %        [~,idx] = findpeaks(wx(f,:));
%        WX(f,idx) = wx(f,idx);
%     end
%     DS = 10;%5 ms bins  
%     WX = movsum(WX,DS,2);
%     WX = WX(:,DS:DS:end);
% %     wx = movmean(wx,DS,2);
%     wx = wx(:,DS:DS:end);
%     tt{ii} = tt{ii}(DS:DS:end);
%     t{ii} = tt{ii};
%     theta_phase = theta_phase(DS:DS:end);
%     wave_phase = wave_phase(DS:DS:end);
%     phase{ii} = phase{ii}(DS:DS:end);
% %     lm = zeros(size(WX,2),25*length(freqVec));
% %     twomax = zeros(size(WX,2),2*length(freqVec));
% %     for f = 1:length(freqVec)
% %         [~,idx] = findpeaks(wx(f,:),'MinPeakProminence',0.025*mean(wx(:)));
% %         for d = 1:size(WX,2)
% %            idx2 = find(,'last'); 
% %         end        
% %         ds = floor(Fs/(2*freqVec(f)));
% %         lwx = lagGen(wx(f,:)',0:ds:26*ds);
% %         lwx = lwx(DS:DS:end,:);
% %         lwx(~islocalmax(lwx,2)) = 0;
% %         lm(:,25*(f-1)+1:25*f) = lwx(:,2:end-1);
% %     end
    
    %save csc data as .csv
    names = cell(1,length(freqVec)+3);
    %     names = cell(1,4*length(freqVec)+1);
    for f = 1:length(freqVec)
        names{f} = strcat('Pow_f',num2str(f));
        %         names{f} = strcat('PICS_f',num2str(f));
        %         names{f+length(freqVec)} = strcat('CICS_f',num2str(f));
        %         names{f+2*length(freqVec)} = strcat('RephasePP_f',num2str(f));
        %         names{f+3*length(freqVec)} = strcat('ImphasePP_f',num2str(f));
    end
    %     names{7*length(freqVec)+1} = 'Time';
    %     names{7*length(freqVec)+2} = 'Position';
    names{length(freqVec)+1} = 'Time';
    names{length(freqVec)+2} = 'Theta_Phase';
    names{length(freqVec)+3} = 'Wave_Phase';
    %     names{4*length(freqVec)+1} = 'Time';
    mt = min(t{ii});
    t{ii} = t{ii} - mt;
    tt{ii} = tt{ii} - mt;
    
    %     rtable = array2table([PICS',CICS',real(phasePP'),imag(phasePP'),tt{ii}'],'VariableNames',names);
    %     rfilename = sprintf('%s%s%s%s',strcat(sessions{ii},'CS\'),'rtable_CSC.csv');
    %     writetable(rtable,rfilename);clear rtable PICS dPICS CICS dCICS phasePP dPP
   
    rtable = array2table([Wx',tt{ii}',theta_phase',wave_phase'],'VariableNames',names);
    rfilename = sprintf('%s%s%s%s',strcat(sessions{ii},'CS\'),'rtable_CSC.csv');
    writetable(rtable,rfilename);clear rtable Wx
    
%     rtable = array2table([WX',tt{ii}',theta_phase',wave_phase'],'VariableNames',names);
%     rfilename = sprintf('%s%s%s%s',strcat(sessions{ii},'CS\'),'rtable_maximan.csv');
%     writetable(rtable,rfilename);clear rtable Wx
    
%     %save lag matrix
%     rtable = array2table(lm);
%     rfilename = sprintf('%s%s%s%s',strcat(sessions{ii},'CS\'),'rtable_lags.csv');
%     writetable(rtable,rfilename);clear rtable Wx
    
%     %save theta asymmetry
%     sym(:,1) = sym(:,1) - mt;
%     names = cell(1,2);
%     names{1} = 'Time';
%     names{2} = 'Asymmetry';
%     rtable = array2table(sym,'VariableNames',names);
%     rfilename = sprintf('%s%s%s%s',strcat(sessions{ii},'CS\'),'sym.csv');
%     writetable(rtable,rfilename);clear theta_stuff
    
    %save tracking data as .csv
    trnames = cell(1,2);
    trnames{1} = 'Time';
    trnames{2} = 'Position';
%     rtable = array2table([t{ii},phase{ii}],'VariableNames',trnames);
    rtable = array2table([tt{ii}',phase{ii}'],'VariableNames',trnames);
    rfilename = sprintf('%s%s%s%s',strcat(sessions{ii},'CS\'),'rtable_tracking.csv');
    writetable(rtable,rfilename);clear rtable
%     %save spiking data as .csv
%     tsnames = cell(1,3);
%         tsnames = cell(1,2*length(freqVec)+2);
%     tsnames{1} = 'CellID';
%     tsnames{2} = 'SpTime';
%     tsnames{3} = 'Theta_Phase';
%     tsnames{4} = 'Wave_Phase';
%     for f = 1:length(freqVec)
%         tsnames{2+f} = strcat('Phase_f',num2str(f));
%         %         tsnames{2+f} = strcat('RephasePP_f',num2str(f));
%         %         tsnames{2+f+length(freqVec)} = strcat('ImphasePP_f',num2str(f));
%     end
%     TS(:,2) = TS(:,2) - mt;
%     TS = [TS PS'];
%     rtable = array2table([TS],'VariableNames',tsnames);
%     rfilename = sprintf('%s%s%s%s',strcat(sessions{ii},'CS\'),'rtable_TS.csv');
%     writetable(rtable,rfilename)
    
    clear CS Wx
%     clear TS
    
end

% %save reward location
% pbins = linspace(-pi,pi,360);
% kappa = 300;
% V = [vel{1};vel{2};vel{3}];
% P = [phase{1};phase{2};phase{3}];
%
% delta = angle(exp(1i*repmat(P,1,length(pbins))).*conj(exp(1i*repmat(pbins,length(P),1))));
% W = exp(kappa*cos(delta));
% D = ones(1,length(V))*exp(kappa*cos(delta));
%
% vmap = V'*W./D;
%
% [~,im] = min(vmap);
% rewloc = pbins(im);
% name = cell(1);
% name{1} = 'rewloc';
%
% dirInfo = dir(parent);
% found = 0;
% for kk=1:size(dirInfo,1)
%     if dirInfo(kk).isdir
%         if strcmp(dirInfo(kk).name,strcat('rewloc','\'))
%             found = 1;
%         end
%     end
% end
% if found==0
%     mkdir(parent,strcat('rewloc','\'));
% end
%
% rtable = array2table(rewloc,'VariableNames',name);
% rfilename = sprintf('%s%s%s%s',parent,strcat('\rewloc','\rewloc'),'.csv');
% writetable(rtable,rfilename)



%%%%%%%%%%%% Other Functions %%%%%%%%%%%%%%

function [dTargets,trackingColour] = decodeTargets(targets)

% Number of samples
numSamp = size(targets,2);

% Allocate memory to the array. 9 fields per sample: X-coord, Y-coord and
% 7 colour flag.
% Colour flag: 3=luminance, 4=rawRed, 5=rawGreen, 6=rawBlue, 7=pureRed,
% 8=pureGreen, 9=pureBlue.
dTargets = int16(zeros(numSamp,50,9));

for ii = 1:numSamp
    for jj = 1:50
        bitField = bitget(targets(jj,ii),1:32);
        if bitField(13)% Raw blue
            % Set the x-coord to the target
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            % Set the y-coord to the target
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,6) = 1;
        end
        if bitField(14) % Raw green
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,5) = 1;
        end
        if bitField(15) % Raw red
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,4) = 1;
        end
        if bitField(16) % Luminance
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,3) = 1;
        end
        if bitField(29) % Pure blue
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,9) = 1;
        end
        if bitField(30) % Puregreen
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,8) = 1;
        end
        if bitField(31) % Pure red
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,7) = 1;
        end
    end
end

% Find out what colours were used in the tracking
trackingColour = zeros(1,7);
if ~isempty(find(dTargets(:,:,3),1)) % Luminance
    trackingColour(1) = 1;
end
if ~isempty(find(dTargets(:,:,7),1)) % Pure Red
    trackingColour(2) = 1;
end
if ~isempty(find(dTargets(:,:,8),1)) % Pure Green
    trackingColour(3) = 1;
end
if ~isempty(find(dTargets(:,:,9),1)) % Pure Blue
    trackingColour(4) = 1;
end
if ~isempty(find(dTargets(:,:,4),1)) % Raw Red
    trackingColour(5) = 1;
end
if ~isempty(find(dTargets(:,:,5),1)) % Raw Green
    trackingColour(6) = 1;
end
if ~isempty(find(dTargets(:,:,6),1)) % Raw Blue
    trackingColour(7) = 1;
end

% Exctracts the individual coordinates for the centre of mass of each
% tracking diode. The red LEDs are assumed to be at the front and the green
% diodes are assumed to be at the back.
function [frontX,frontY,backX,backY] = extractPosition(targets,tracking)

ind = find(tracking(2:end));
if length(ind) <= 1
    % Need at least two colours to get head direction
    disp('ERROR: To few LED colours have been tracked. Not possible to find head direction')
    frontX = NaN;
    frontY = NaN;
    backX = NaN;
    backY = NaN;
    return
else
    if ~tracking(2) && ~tracking(5)
        disp('ERROR: Red LED has not been tracked')
        frontX = NaN;
        frontY = NaN;
        backX = NaN;
        backY = NaN;
        return
    end
    if ~tracking(3) && ~tracking(6)
        disp('ERROR: Green LED has not been tracked')
        frontX = NaN;
        frontY = NaN;
        backX = NaN;
        backY = NaN;
        return
    end
end

% Number of samples in the data
numSamp = size(targets,1);

% Allocate memory for the arrays
frontX = zeros(1,numSamp);
frontY = zeros(1,numSamp);
backX = zeros(1,numSamp);
backY = zeros(1,numSamp);

% Exctract the front coordinates (red LED)
if tracking(2) && ~tracking(5)
    % Pure red but not raw red
    for ii = 1:numSamp
        ind = find(targets(ii,:,7));
        if ~isempty(ind)
            frontX(ii) = mean(targets(ii,ind,1));
            frontY(ii) = mean(targets(ii,ind,2));
        end
    end
end
if ~tracking(2) && tracking(5)
    % Not pure red but raw red
    for ii = 1:numSamp
        ind = find(targets(ii,:,4));
        if ~isempty(ind)
            frontX(ii) = mean(targets(ii,ind,1));
            frontY(ii) = mean(targets(ii,ind,2));
        end
    end
end
if tracking(2) && tracking(5)
    % Both pure red and raw red
    for ii = 1:numSamp
        ind = find(targets(ii,:,7) | targets(ii,:,4));
        if ~isempty(ind)
            frontX(ii) = mean(targets(ii,ind,1));
            frontY(ii) = mean(targets(ii,ind,2));
        end
    end
end

% Exctract the back coordinates (green LED)
if tracking(3) && ~tracking(6)
    % Pure green but not raw green
    for ii = 1:numSamp
        ind = find(targets(ii,:,8));
        if ~isempty(ind)
            backX(ii) = mean(targets(ii,ind,1));
            backY(ii) = mean(targets(ii,ind,2));
        end
    end
end
if ~tracking(3) && tracking(6)
    % Not pure green but raw green
    for ii = 1:numSamp
        ind = find(targets(ii,:,5));
        if ~isempty(ind)
            backX(ii) = mean(targets(ii,ind,1));
            backY(ii) = mean(targets(ii,ind,2));
        end
    end
end
if tracking(3) && tracking(6)
    % Both pure green and raw green
    for ii = 1:numSamp
        ind = find(targets(ii,:,8) | targets(ii,:,5));
        if ~isempty(ind)
            backX(ii) = mean(targets(ii,ind,1));
            backY(ii) = mean(targets(ii,ind,2));
        end
    end
end

% Estimates lacking position samples using linear interpolation. When more
% than timeTreshold sec of data is missing in a row the data is left as
% missing.
%
% Raymond Skjerpeng 2006.
function [theta,r] = interporPos(x,y,timeTreshold,sampRate)

theta = unwrap(atan2(y,x));
r = sqrt(x.^2+y.^2);

% Turn off warning
% warning('off','MATLAB:divideByZero');

% Number of samples that corresponds to the time threshold.
sampTreshold = floor(timeTreshold * sampRate);

% number of samples
numSamp = length(x);
% Find the indexes to the missing samples
% temp1 = 1./x;
% indt1 = isinf(temp1);
ind = isnan(x);
ind2 = find(ind==1);
% Number of missing samples
N = length(ind2);

if N == 0
    % No samples missing, and we return
    return
end

change = 0;

% Remove NaN in the start of the path
if ind2(1) == 1
    change = 1;
    count = 0;
    while 1
        count = count + 1;
        if ind(count)==0
            break
        end
    end
    theta(1:count) = theta(count);
    r(1:count) = r(count);
end

% Remove NaN in the end of the path
if ind2(end) == numSamp
    change = 1;
    count = length(x);
    while 1
        count = count - 1;
        if ind(count)==0
            break
        end
    end
    theta(count:numSamp) = theta(count);
    r(count:numSamp) = r(count);
end

if change
    % Recalculate the missing samples
    %     temp1 = 1./r;
    %     indt1 = isinf(temp1);
    ind = isnan(r);
    % Missing samples are where both x and y are equal to zero
    ind2 = find(ind==1);
    % Number of samples missing
    N = length(ind2);
end

for ii = 1:N
    % Start of missing segment (may consist of only one sample)
    start = ind2(ii);
    % Find the number of samples missing in a row
    count = 0;
    while 1
        count = count + 1;
        if ind(start+count)==0
            break
        end
    end
    % Index to the next good sample
    stop = start+count;
    if start == stop
        % Only one sample missing. Setting it to the last known good
        % sample
        theta(start) = theta(start-1);
        r(start) = r(start-1);
    else
        if count < sampTreshold
            % Last good position before lack of tracking
            theta1 = theta(start-1);
            r1 = r(start-1);
            % Next good position after lack of tracking
            theta2 = theta(stop);
            r2 = r(stop);
            % Calculate the interpolated positions
            theta(start:stop) = interp1([1,2],[theta1,theta2],1:1/count:2);
            r(start:stop) = interp1([1,2],[r1,r2],1:1/count:2);
            % Increment the counter (avoid estimating allready estimated
            % samples)
            ii = ii+count;
        else
            % To many samples missing in a row and they are left as missing
            ii = ii+count;
        end
    end
end
theta = wraptopi(theta);

% Calculates the direction of the head stage from the two set of
% coordinates. If one or both coordinate sets are missing for one samle the
% direction is set to NaN for that sample. Direction is also set to NaN for
% samples where the two coordinate set are identical. Returns the
% direction in degrees
function direct = headDirection(frontX,frontY,backX,backY)

% Number of position samples in data set
N = length(frontX);
direct = zeros(N,1);

for ii = 1:N
    
    if frontX(ii)==0 || backX(ii)==0 || isnan(frontX(ii)) || isnan(backX(ii))
        % One or both coordinates are missing. No angle.
        direct(ii) = NaN;
        continue
    end
    
    % Calculate the difference between the coordinates
    xd = frontX(ii) - backX(ii);
    yd = frontY(ii) - backY(ii);
    
    if xd==0
        if yd==0
            % The two coordinates are at the same place and it is not
            % possible to calculate the angle
            direct(ii) = NaN;
            continue
        elseif yd>0
            direct(ii) = 90;
            continue
        else
            direct(ii) = 270;
            continue
        end
    end
    if yd==0
        if xd>0
            % Angle is zero
            continue
        else
            direct(ii) = 180;
            continue
        end
    end
    
    if frontX(ii)>backX(ii) && frontY(ii)>backY(ii)
        % Angle between 0 and 90 degrees
        direct(ii) = atan(yd/xd) * 360/(2*pi);
        
    elseif frontX(ii)<backX(ii) && frontY(ii)>backY(ii)
        % Angle between 90 and 180 degrees
        direct(ii) = 180 - atan(yd/abs(xd)) * 360/(2*pi);
        
    elseif frontX(ii)<backX(ii) && frontY(ii)<backY(ii)
        % Angle between 180 and 270 degrees
        direct(ii) = 180 + atan(abs(yd)/abs(xd)) * 360/(2*pi);
        
    else
        % Angle between 270 and 360 degrees
        direct(ii) = 360 - atan(abs(yd)/xd) * 360/(2*pi);
    end
end

function [sptimes] = spikeTimes(ts,loct)

N = length(ts);
sptimes = zeros(N,4);

for ii = 1:N
    tdiff = loct-ts(ii);
    sptimes(ii,1) = find(tdiff<=0,1,'last');
    sptimes(ii,2) = tdiff(sptimes(ii,1));
    sptimes(ii,3) = find(tdiff>0,1,'first');
    sptimes(ii,4) = tdiff(sptimes(ii,3));
end

