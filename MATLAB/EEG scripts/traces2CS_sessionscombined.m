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

function [CS] = traces2CS_sessionscombined(inFile,freqVec,width,numperms,pdim)

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

pdim = numchannels;
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
nbands = sum(freqVec>=6&freqVec<=12); %bands used for theta pca
thbandsf = [];
N = zeros(numsessions,1);
MCf = zeros(length(freqVec),numchannels,numchannels);
t{ii} = cell(numsessions,1);
phase{ii} = cell(numsessions,1);
vel{ii} = cell(numsessions,1);

for ii = 1:numsessions
     disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));
    
    % Check if subdir for storing images are present. If not, it is
    % created
    dirInfo = dir(sessions{ii});
    found = 0;
    for kk=1:size(dirInfo,1)
        if dirInfo(kk).isdir
            if strcmp(dirInfo(kk).name,strcat('CS_comb','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(sessions{ii},strcat('CS_comb','\'));
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
    frontX(frontX==0&frontY==0) = NaN;frontY(frontX==0&frontY==0) = NaN;
    backX(backX==0&backY==0) = NaN;backY(backX==0&backY==0) = NaN;
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
    phase{ii} = circ_mean([fphase;bphase]);
    %smooth and get derivatives
    phase{ii} = medfilt1(unwrap(phase{ii}),5);
    [phase{ii},~] = localfit_ad(t{ii},phase{ii},2,0.04);
    phase{ii} = medfilt1(unwrap(phase{ii}),5);
    [~,vel{ii}] = localfit_ad(t{ii},phase{ii},2,0.175);
    t{ii}(isnan(vel{ii})) = [];phase{ii}(isnan(vel{ii})) = [];vel{ii}(isnan(vel{ii})) = [];
    vel{ii} = medfilt1(vel{ii},5);
    phase{ii} = wraptopi(phase{ii});
    phase{ii} = reshape(phase{ii},length(phase{ii}),1);
    vel{ii} = reshape(vel{ii},length(vel{ii}),1);
    vel{ii} = vel{ii}*D/2;
    
    %load wavelet transforms
    Wx = zeros(length(freqVec),length(t{ii}),numchannels);    
    for jj = 1:numchannels
        % Load data from the .ncs files, make plots, and store them
        xfile = [sessions{ii},channels{jj}];
        [samples,~,tt, Fs, bv, ~] = loadEEG2(xfile);
        x = bv*samples;
        
        wx = traces2Wx(x,freqVec,Fs,width);
        if jj == 1  
            thbandst = zeros(length(tt),numchannels*nbands);
            ttidx = [];
            for it = 1:length(t{ii})
                [~,idx] = min((tt-t{ii}(it)).^2);
                ttidx = [ttidx idx];
            end            
        end

        Wx(:,:,jj) = wx(:,ttidx);
        thbandst(:,(jj-1)*nbands+1:jj*nbands) = wx(freqVec>=6&freqVec<=12,:).';
        clear wx
    end   

    thbandsf = [thbandsf;thbandst]; clear thbandst
    
    %compute cross spectra
    CS = zeros(length(freqVec),length(t{ii}),numchannels,numchannels);
    for jj = 1:numchannels
        CS(:,:,:,jj) = bsxfun(@times,Wx,conj(Wx(:,:,jj)));
    end
    
    MCf = MCf+squeeze(sum(CS,2));
    N(ii) = size(CS,2); clear CS
end

MCf = MCf/sum(N);
MCf = permute(MCf,[2 3 1]);
thc = pca(thbandsf,'Centered',false);clear thbandsf

for ii = 1:numsessions   

    %load wavelet transforms and do PCA on theta
%     Wx = cell(numchannels,1);
    Wx = zeros(length(freqVec),length(t{ii}),numchannels);
    nbands = sum(freqVec>=6&freqVec<=12);    
    for jj = 1:numchannels
        % Load data from the .ncs files, make plots, and store them
        xfile = [sessions{ii},channels{jj}];
        [samples,~,tt, Fs, bv, ~] = loadEEG2(xfile);
        x = bv*samples;
        
        wx = traces2Wx(x,freqVec,Fs,width);
        if jj == 1  
            thbands = zeros(length(tt),numchannels*nbands);
            ttidx = [];
            for it = 1:length(t{ii})
                [~,idx] = min((tt-t{ii}(it)).^2);
                ttidx = [ttidx idx];
            end            
        end

        Wx(:,:,jj) = wx(:,ttidx);
        thbands(:,(jj-1)*nbands+1:jj*nbands) = wx(freqVec>=6&freqVec<=12,:).';
        clear wx
    end   

    theta = sum(repmat(thc(:,1).',length(t{ii}),1).*thbands(ttidx,:),2); clear thbands
    thetaphase = angle(theta);clear theta
    
    %compute cross spectra
    CS = zeros(length(freqVec),length(t{ii}),numchannels,numchannels);
    for jj = 1:numchannels
        CS(:,:,:,jj) = bsxfun(@times,Wx,conj(Wx(:,:,jj)));
    end
    
%     start = round(linspace(310,size(Wx,2)-310,numperms));
%     pMC = zeros(length(freqVec),numchannels,numchannels,numperms); 
%     pcs = zeros(length(freqVec),length(t{ii}),numchannels,numchannels);
%     for p = 1:numperms
%         for jj = 1:numchannels
%             pcs(:,:,:,jj) = bsxfun(@times,Wx,conj([Wx(:,start(p):end,jj) Wx(:,1:start(p)-1,jj)]));            
%         end
%         pMC(:,:,:,p) = squeeze(sum(pcs,2));
%     end
%     clear Wx pcs
%     for jj = 1:numchannels
%        pMC(:,jj,jj,:) = repmat(MC(:,jj,jj),1,1,1,numperms); 
%     end
    
%     pMC = permute(pMC,[2 3 1 4]);
    GIM = zeros(length(freqVec),1);
%     pGIM = zeros(length(freqVec),numperms);
    PCS = zeros(length(freqVec),length(t{ii}));
%     GIMn = zeros(numchannels,length(freqVec));
    alpha = zeros(pdim,length(freqVec)); beta = alpha;
    re = zeros(pdim,pdim,length(freqVec));
%     reg = 0.001;
    for f = 1:length(freqVec)
%         keyboard
        mc = MCf(:,:,f) - tril(MCf(:,:,f),-1) + tril(MCf(:,:,f)',-1);
%         scaler = mean(mean(abs(mc)));
%         mc = mc./scaler; %avoid small number numerical issues
%         cit=imag(mc)*imag(mc)'+eye(numchannels,numchannels)*mean(diag(imag(mc)*imag(mc)'))*10^(-10);  % Regularization
%         [u,~,~] = svd(cit);
%         mc=u(:,1:pdim)'*mc*u(:,1:pdim); % Dimensionality reduction
%         mc=mc+eye(pdim,pdim)*(mean(diag(real(mc)))*reg); % Regularization
        re(:,:,f) = real(mc)^-0.5;
        rho = re(:,:,f)*imag(mc)*re(:,:,f);
        [V,~] = eig(rho*rho');
        alpha(:,f) = V(:,end);
        beta(:,f) = V(:,end-1);
%         CS(f,:,:,:) = CS(f,:,:,:)/scaler;
        for it=1:length(t{ii})
%             pdcs = u(:,1:pdim)'*squeeze(CS(f,it,:,:))*u(:,1:pdim);
%             pdcs = pdcs+eye(pdim,pdim)*(mean(diag(real(mc)))*reg); % Regularization
            PCS(f,it) = alpha(:,f)'*re(:,:,f)*squeeze(CS(f,it,:,:))*re(:,:,f)*beta(:,f);
        end   
%         GIMn(:,f) = eig(real(mc)\imag(mc)*inv(real(mc))*imag(mc)');        
        GIM(f) = 0.5*trace(real(mc)\imag(mc)*inv(real(mc))*imag(mc)');
%         for p = 1:numperms
%             mc = pMC(:,:,f,p) - tril(pMC(:,:,f,p),-1) + tril(pMC(:,:,f,p)',-1);
%             mc = mc./mean(mean(abs(mc)));
%             pGIM(f,p) = 0.5*trace(real(mc)\imag(mc)*inv(real(mc))*imag(mc)');
%         end
    end
    close all
%     keyboard
    subplot(2,2,1)
    ph1 = plot(freqVec,GIM,'-k','LineWidth',2);hold on
    ph2 = plot(freqVec,abs(mean(imag(PCS),2)).^2,'b','LineWidth',2);
%     ph3 = plot(freqVec,mean(pGIM,2),'-r','LineWidth',2); axis square
    %hleg1 = legend('GIM','1st Component','Shuffled');
    sh=subplot(2,2,2);
    p=get(sh,'position');
    lh=legend(sh,[ph1;ph2],'GIM','1st Component');
    set(lh,'position',p);
    axis(sh,'off');
    %set(hleg1,'Location','northoutside','Orientation','horizontal')
    subplot(2,2,1); hold on
    xlabel('Frequency (Hz)');ylabel('GIM');
    ylim([0 max(max([GIM abs(mean(imag(PCS),2)).^2]))])
%     [ci] = prctile(pGIM',[2.5 97.5]);
%     shadedplot(freqVec, ci(2,:), ci(1,:),[1 0 0],[1 0.25 0.25]);
    
    subplot(2,2,3)
    plot(freqVec,abs(mean(imag(PCS),2)).^2./GIM,'-k','LineWidth',2); axis square
    xlabel('Frequency (Hz)')
    ylabel('Proportion of GIM')
%     zGIM = (GIM - mean(pGIM,2))./std(pGIM,0,2);
%     plot(freqVec,zGIM,'-k','LineWidth',2); axis square
%     xlabel('frequency (Hz)');ylabel('GIM (z-score from shuffled)');
%     ylim([min(zGIM) max(zGIM)])
    
    figImage = sprintf('%s%s%s%s',sessions{ii},strcat('CS_comb','\GIM'),'.fig');
    bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('CS_comb','\GIM'),'.bmp');
    saveas(gcf,figImage,'fig');
    saveas(gcf,bmpImage,'bmp');
    
    %save data as struct
    data.freqVec = freqVec; data.CS = CS; data.PCS = PCS; data.alpha = alpha; data.time = t{ii}; 
    data.beta = beta; data.phase = phase{ii}; data.vel = vel{ii}; data.thetaphase = thetaphase; 
    filename = sprintf('%s%s%s%s',sessions{ii},strcat('CS_comb','\CS'),'.mat');
    save(filename,'-struct','data');
    clear CS data
    
    data.freqVec = freqVec; data.GIM = GIM; %data.pGIM = pGIM; %data.zGIM = zGIM;
    filename = sprintf('%s%s%s%s',sessions{ii},strcat('CS_comb','\data'),'.mat');
    save(filename,'-struct','data');
    clear data
end

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


