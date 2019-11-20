function [modes] = GMD(inFile,freqVec,pbins)

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

V = [];

for ii = 1:numsessions
    disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));
    % Check if subdir for storing images are present. If not, it is
    % created
    dirInfo = dir(sessions{ii});
    found = 0;
    for kk=1:size(dirInfo,1)
        if dirInfo(kk).isdir
            if strcmp(dirInfo(kk).name,strcat('GMD','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(cd,strcat('GMD','\'));
    end
    
    % Get position data
    fieldSelection(1) = 1; % Timestamps
    fieldSelection(2) = 1; % Extracted X
    fieldSelection(3) = 1; % Extracted Y
    fieldSelection(4) = 0; % Extracted Angel
    fieldSelection(5) = 0; % Targets
    fieldSelection(6) = 0; % Points
    % Do we return header 1 = Yes, 0 = No.
    extractHeader = 0;
    % 5 different extraction modes, see help file for Nlx2MatVt
    extractMode = 1; % Extract all data
    
    %compute running speed
    % Get position data
    file = strcat(sessions{ii},'vt1.nvt');
    [t{ii},x,y] = Nlx2MatVT(file,fieldSelection,extractHeader,extractMode);
    % Convert timestamps to seconds
    t{ii} = t{ii}'/1000000;
    %     fidx = x==0&y==0;
    
    %interpolate positions
    x = pchip(t{ii}(x~=0),x(x~=0),t{ii});
    y = pchip(t{ii}(y~=0),y(y~=0),t{ii});
    
    %center track
    %     nb = 250;
    %     xymap = zeros(nb,nb);
    %     v = 0.1;
    %     mx = linspace(min(x),max(x),nb);
    %     my = linspace(min(y),max(y),nb);
    %     for xx = 1:length(mx)
    %         dx = mx(xx)-x;
    %         for yy = 1:length(my)
    %             dy = my(yy)-y;
    %             xymap(xx,yy) = nansum(v^2/sqrt(2*pi)*exp(-0.5*(dx.*dx*v^2+dy.*dy*v^2)));
    %         end
    %     end
    %     thmap = xymap+0.1;thmap(xymap>quantile(xymap(:),0.65))=0;
    
    %     regmax = imregionalmax(xymap,4);
    %     [xx,yy] = find(regmax);
    %     [z, a, b, alpha] = fitellipse([mx(xx);my(yy)]); %'linear' to speed up
    %     %Translate
    %     x = x-z(1);
    %     y = y-z(2);
    %     %Rotate
    %     Q = [cos(-alpha), -sin(-alpha); sin(-alpha) cos(-alpha)];
    %     temp = Q*[x';y'];
    %     %Scale
    %     D = 100;
    %     temp(1,:) = 0.5*D/a*temp(1,:);
    %     temp(2,:) = 0.5*D/b*temp(2,:);
    %     %Rotate back to orginal orientation
    %     Q = [cos(alpha), -sin(alpha); sin(alpha) cos(alpha)];
    %     temp = Q*temp;
    %     x = temp(1,:);y = temp(2,:);
    
    [xc yc R] = circfit(x,y);
    x = x - xc;
    y = y - yc;
    badidx = find(diff(t{ii})<0)+1;
    x(badidx)=[];y(badidx)=[];
    t{ii}(badidx)=[];
    
    %convert position to phase, smooth, and get velocity
    D = 100;
    pos = -unwrap(atan2(y,x));
    span = 4/(max(t{ii})-min(t{ii}));
    pos = smooth(pos,span,'lowess');
    vel = abs(diff(pos))./diff(t{ii});
    vel = [vel(1);vel];
    span = 1/(max(t{ii})-min(t{ii}));
    vel = exp(smooth(real(log(vel)),span,'lowess'))*D/2;
    V = [V;vel];
%     t{ii}(vel<10) = [];
end

% Load data from the .ncs files, make plots, and store them
for jj=1%:numchannels
    TFR = [];TP = [];X = [];
    for ii = 1:numsessions
        xfile = [sessions{ii},channels{jj}];
        [x,~,tt, Fs, bv, ir] = loadEEG2(xfile);
        x = bv*x;keyboard
%         xx = x(10*2000:60*2000)';
        
        %get high/bandpass
        HP = ffthighpass(x,Fs,5,7);
        BP = fftlowpass(HP,Fs,11,13);
        
        %get monotonic phase
        [ppks,pidx] = findpeaks(-BP,'MinPeakDistance',Fs/13);
        phase = pchip(pidx,1:length(pidx),1:length(BP))';
        
        %get instantateous ampliude
        amp = abs(hilbert(BP));
        
        %construct matrix of theta cycles
        ngrid = 350;
        pgrid = 0:1/ngrid:1-1/ngrid;   
        HPint = pchip(phase,HP./amp,1:1/ngrid:floor(max(phase))-1/ngrid)';        
        HPint = reshape(HPint,ngrid,length(HPint)/ngrid);
        Velint = pchip(phase,vel',1:1/ngrid:floor(max(phase))-1/ngrid)';        
        Velint = reshape(Velint,ngrid,length(Velint)/ngrid);
        mindev = max(abs(Velint));
        HPint = HPint(:,mindev>5);
        
        %  -------------  set up fourier basis  ---------------------------
        nbasis     = 150;
        prange = [min(pgrid) max(pgrid)];
        fbasis = create_fourier_basis(prange, nbasis);
        
        %  ----  set up the harmonic acceleration operator  -------
        Lbasis  = create_constant_basis(prange);  %  create a constant basis
        Lcoef   = [0,(2*pi/1)^2,0];    %  set up three coefficients
        wfd     = fd(Lcoef,Lbasis);      % define an FD object for weight functions
        wfdcell = fd2cell(wfd);          % convert the FD object to a cell object
        harmaccelLfd = Lfd(3, wfdcell);  %  define the operator object
        
        % set cross validation parameters
        numtest = round(0.1*size(HPint,1));
        testidx = randsample(size(HPint,1),numtest);
        HPtrain = HPint;
        HPtrain(:,testidx) = [];
        n = size(HPtrain,2);
        lambda   = 10.^-(9:-0.3:-1);        
        maxharm = 20;
        CV = zeros(length(lambda),maxharm);
        
        % perform cross validation on lambda and nharm
        for l = 1:length(lambda)
            %setup smoothing objects
            fdParobj = fdPar(fbasis, harmaccelLfd, lambda(l));
            [thetafd, df, ~, ~, ~] = smooth_basis(pgrid, HPtrain, fdParobj);
            fdParobj = fdPar(thetafd, harmaccelLfd, lambda(l)); 
            %do pca
            thetapcastr = pca_fd(thetafd, maxharm, fdParobj);
            %get mean and harmonic functions
            meanmat = squeeze(eval_fd(thetapcastr.meanfd, pgrid));
            harmmat = eval_fd(thetapcastr.harmfd, pgrid);            
            %loop through number of harmonic components
            fitm = repmat(meanmat,1,numtest);
            for m = 1:maxharm  
                SSEc = zeros(1,numtest);
                %loop through test cycles
                for c = 1:numtest
%                     score = sum((HPint(:,[testidx(c)-1:testidx(c)-1 testidx(c)+1:testidx(c)+1])-meanmat)'*harmmat(:,m))/(ngrid*4); 
                    fitm(:,c) = fitm(:,c);% + score*harmmat(:,m); 
                    SSEc(c) = sum((HPint(:,testidx(c)) - fitm(:,c)).^2);
                end
                CV(l,m) = n/(n-df)*sum(SSEc)/(numtest*ngrid)/(n-df);
            end
        end
        plot(log10(lambda),CV(:,1),'ok')
            
        lambda = 10^-(6);
        nharm = 10;
        %setup smoothing objects
        fdParobj = fdPar(fbasis, harmaccelLfd, lambda);
        [thetafd, ~, ~, ~, ~] = smooth_basis(pgrid, HPint, fdParobj);
        fdParobj = fdPar(thetafd, harmaccelLfd, lambda);
        %do pca 
        thetapcastr = pca_fd(thetafd, nharm, fdParobj,0);
        plot_pca_fd(thetapcastr, 1)
        
        %get harmonic functions
        phaseint = pgrid(round(ngrid*mod(phase(1),1))):1/ngrid:(pgrid(round(ngrid*mod(phase(end),1)))+floor(phase(end)));
        ampint = pchip(phase,amp,phaseint)';
        HPrec = pchip(phase,HP./amp,phaseint)';
        BPrec = pchip(phase,BP,phaseint)';
        for p = 1:length(phaseint)
           [~,idx] = min(abs(pgrid-mod(phaseint(p),1))); 
           phaseint(p) = pgrid(idx);
        end
        harmmat = eval_fd(thetapcastr.harmfd, phaseint);
        meanmat = squeeze(eval_fd(thetapcastr.meanfd, phaseint));
        
        %reconstruct signal
        winsize = 1*ngrid;
        window = ones(winsize,1)/(winsize);
        innermat = harmmat.*(HPrec-meanmat);
        coefmat = zeros(size(innermat));
        for c = 1:nharm
           tmpc = convfft(innermat(:,c),window);
           coefmat(:,c) = tmpc(ceil(winsize/2):length(tmpc)-floor(winsize/2));
        end
        fit = meanmat + sum(harmmat(:,1:10).*coefmat(:,1:10),2);
        
        plot(HPrec.*ampint);hold on
        plot(BPrec,'-b','LineWidth',2);plot(fit.*ampint,'-r','LineWidth',2);
        xlim([1070000 1100000])
        
%         %set sswp parameters
%         N = numel(xx);
%         eps = 1e-2;
%         res = 0.2*N/Fs;
%         freq_range = [0 100*N/Fs];
%         NG = N;
%         t_sc = 0.85;
%         xt = 0:1/N:(1-1/N);
%         red = 10;
%         rad = 1;
%         %sswp
%         [T_f, coef, kk] = ss_wp1_fwd(xx,1,1,1,xt,NG,freq_range(2),freq_range(1),rad,1,t_sc,red,eps,res,0);
%         %determine component frequency range
%         sp = smooth(log10(linspace(0,150,600)),nanmean(T_f,2),0.05,'lowess');
%         [~,pidx] = findpeaks(sp);
%         [~,vidx] = findpeaks(-sp);
%         edidx = vidx(find(vidx>pidx(1),1,'first'));
%         stidx = 2*pidx(1)-edidx;
%         %extract component in given ranges
%         thre = 0;
%         C = 1;
%         max_num = 1;
%         pct = 0.01;
%         T_temp = cell(1,1);
%         T_temp{1} = zeros(size(T_f));
%         T_temp{1}(stidx:edidx,:) = T_f(stidx:edidx,:);        
%         [cluster, freq, lb, ub] = freq_selection(T_temp{1}, 1, eps*10, C, max_num, thre, res, pct, freq_range);
%         [mode, ~] = ss_wp1_invT(cluster, coef, kk, 1, N, freq_range(2), freq_range(1), rad, 1, t_sc, res,0);        
%         %adaptive high pass
%         lowclust = cell(1,1);
%         lowclust{1} = zeros(size(T_f));
%         idx = 1:size(T_f,1);
%         for i = 1:size(T_f,2)
%             lowclust{1}(idx<lb(i)/res,i) = T_f(idx<lb(i)/res,i);
%         end
%         [lowmode, ~] = ss_wp1_invT(lowclust, coef, kk, 1, N, freq_range(2), freq_range(1), rad, 1, t_sc, res,0);
%         hp = x - lowmode;
%         %estimate instantaneous info
%         freq = filterFreq(freq,1);
%         lb = filterFreq(lb,1);
%         ub = filterFreq(ub,1);
%         ins_amplt = amplt_est(1, mode);
%         ins_pre_phase = pre_phase_est(freq,1/N);
%         %correct phases
%         peaks = peakDetection(mode,freq);
%         ins_pre_phase = phaseShift(ins_pre_phase,peaks);

        %save table for use in R
        names{1} = 'x';
        names{2} = 'phase';
        names{3} = 'time';
        names{4} = 'vel';
        names{5} = 'amp';

        rtable = array2table([((HP-BP)./amp) mod(phase,1) (1:length(x))' vel' amp],'VariableNames',names);
        rfilename = sprintf('%s%s%s%s','C:\Users\Brian\Documents\R\inputR.csv');
        writetable(rtable,rfilename)
        %run mgcv R script via system command
        !R CMD BATCH C:\Users\Brian\Documents\R\mgcvGMD.r
        %load output from R
        rtable = csvread(sprintf('%s%s%s%s','C:\Users\Brian\Documents\R\outputR.csv'),1,0);
        [newT_f, newcoef, newkk] = ss_wp1_fwd(hp-rtable.*ins_amplt,1,1,1,xt,NG,freq_range(2),freq_range(1),rad,1,t_sc,red,eps,res,0);
    end
    disp('Make plots and store them to files');
    disp(sprintf('%s%i',' CSC ',jj,' of ',numchannels));
    file = [sessions{ii},channels{jj}];
    
    kappa = 100;
    delta = angle(exp(1i*repmat(newTP',1,length(pbins))).*conj(exp(1i*repmat(pbins,length(TP),1))));
    W = exp(kappa*cos(delta));
    D = ones(1,length(TP))*exp(kappa*cos(delta));
%     pfr_mean = TFR*W./D;
    bcycle = ((HP'-new))*W./D;
    pp = csape(pbins,bcycle,'periodic');
    new = new + fnval(pp,newTP);
    DTAS = hilbert(new);
    newTP = angle(DTAS);   
    

    figure(1);
    col = max(abs([min(min(pfr_mean)) max(max(pfr_mean))]));
    surf(pbins,freqVec,zeros(size(pfr_mean)),pfr_mean);colorbar;view([0 90]);axis tight;axis square;shading interp;
    caxis([-col col]);
    set(gca,'Ytick',0:10:100,'yticklabel',0:10:100,'YScale','log');
    %         imagesc(phasebins,freqVec,pfr_mean); axis xy;colorbar;
    xlabel('phase (rads)');ylabel('frequency (Hz)');
    title(strcat('Phase Frequency Representation (z-scored): ',channels{jj}(1:end-4)));
    hold on
    ncycle = cycle-min(cycle);ncycle=ncycle/max(ncycle);
    plot(pbins,80*ncycle+min(freqVec),'k');hold off

    bmpImage = sprintf('%s%s%s%s',strcat('PFR','\'),channels{jj}(1:end-4),'.bmp');
%     saveas(gcf,bmpImage,'bmp');
end

