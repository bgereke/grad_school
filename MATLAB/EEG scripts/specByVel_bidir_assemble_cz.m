function [vel,subvel,sample,vr,spechist_left,spechist_right,F,runData,S_abs,Sz_abs,S_rlt,Sz_rlt,S_sum]...
    =specByVel_bidir_assemble_cz(infile, odir, csclist,badlap,limit,sigbin,plot_fs,plotcolor)
%%----------------------------------------------------
%	pecByVel_bidir_assemble_cz -- combine the laps across all the sessions
%		inputs:
%			infile:  the name and path of the standard infile.txt
%			odir: the name of the base directory for file outputs
%			csclist: a vector listing the number portion of the csc files to look at (eg. [1 3 12] to look at csc1, csc3 and csc12
%           badlap: ndirs*1 cells, each cell is the number vector of the bad laps
%           limit: to find the run-by-run laps with limit at ends
%           sigbin: the sign of speed bin, sigbin=0 limear speed bin 0:5:100 cm/s, sigbin=1 log2 speed bin
%           plot_fs: plot_fs=[minf,maxf],plot from the min to the max frequency
%           plotsolor: colorbar range
%
%		outputs:
%			vel: cell array (directory) of the raw velocity vector
%			subvel: cell array (directory) of the smoothed velocity vector
%			gamma_subsample: cell array (directory, csc count) of high-pass filtered EEG traces, 
%			Vr: the velocity bins
%			spechist: cell array (directory, csc count) of the normalized, binned spectrogram histograms
%%----------------------------------------------------

%% Pre-process and load files
[oudir, nlines, ndirs, ttfile, ncells, dataloc] = dataPrep(infile, odir);
global inp;
if isempty(inp)
	inp = input('Do you want to interpolate missing position data? (Y/N): ','s');	%check that the user wishes to continue with interpolation
	while( strcmp(lower(inp),'y')~=1 && strcmp(lower(inp), 'n') ~=1)  % lower('str') returnsthe string formed by converting any uppercase characters in str tothe corresponding lowercase characters and leaving all other charactersunchanged. 
		inp = input('Do you want to continue with interpolation? (Y/N): ','s');
	end
end
global dir_delim;
if ispc
	dir_delim='\';
else
	dir_delim = '/';
end

% sigbin=0 limear speed bin 0:5:100 cm/s, sigbin=1 log2 speed bin
if sigbin==0
    % linear vr
    speedbin=5; % cm/s
    nbin=100/speedbin;
    vr=linspace(0,100,nbin+1); %velocity bins in cm/s   % linspace(a,b,n) generatesa row vector y of n points linearlyspaced between and including a and b.For n < 2, linspace returns b.
end
if sigbin==1
    % log2 vr
    vr=2.^[0:0.35:7]*0.75;
    vr=[0,vr];
    nbin=length(vr);
end

neeg=length(csclist);
spechist_left=cell(1,neeg);
spechist_left_z=cell(1,neeg);
bincount_left=cell(1,neeg);
spechist_right=cell(1,neeg);
spechist_right_z=cell(1,neeg);
bincount_right=cell(1,neeg);
spechist_both=cell(1,neeg);
spechist_both_z=cell(1,neeg);
bincount_both=cell(1,neeg);

S_abs=cell(ndirs,neeg);
Sz_abs=cell(ndirs,neeg);
S_rlt=cell(ndirs,neeg);
Sz_rlt=cell(ndirs,neeg);

padval=2;  %additional padding orders of two in the spectrogram

if strcmp(lower(inp),'y')
    outdir = strcat(oudir,dir_delim,'interp',dir_delim);
else
    outdir = strcat(oudir,dir_delim,'nointerp',dir_delim);
end
dirnam=outdir;
if ~isdir(dirnam)
    mkdir(dirnam);
end
subdir1=strcat(dirnam,dir_delim,'Figures',dir_delim);
if ~isdir(subdir1)
    mkdir(subdir1);
end
subdir2=strcat(dirnam,dir_delim,'Assembles',dir_delim);
if ~isdir(subdir2)
    mkdir(subdir2);
end

runData=cell(ndirs,1);

%% get sessions
for aa=1:ndirs  % ndirs is the number of sessions
   % get velocity vector, bin into [velwin] ms chunks
    scale=0.26;
    vfs=29.97;  % ======= the sample frequency of video ========
    timelimit=120;   % ======= set the highest speed, unit is cm/s ========
	[t,x,y] = loadPos_rescaled_cz(strcat(dataloc{aa},dir_delim,'VT1.nvt'),scale,scale,vfs);
    y = max(y)-y;
    
    % Speed
	vel{aa}=speed2D(x,y,t); %velocity in cm/s
	vel{aa}(vel{aa}>=timelimit) = 0.5*(vel{aa}(circshift((vel{aa}>=timelimit),-3)) + vel{aa}(circshift((vel{aa}>=timelimit),3)));
    % B = circshift(A,shiftsize) circularlyshifts the values in the array, A, by shiftsize elements.
    % shiftsize isa vector of integer scalars where the n-th elementspecifies the shift amount for the n-th dimensionof array A. 
    % If an element in shiftsize ispositive, the values of A are shifted down (orto the right). 
    % If it is negative, the values of A areshifted up (or to the left). If it is 0, the valuesin that dimension are not shifted.
	nWin=floor(length(vel{aa})/6);
    subvel{aa}=zeros(1,nWin);  % average every 6 positions data
	for ii=1:nWin
		subvel{aa}(ii)=mean(vel{aa}((ii-1)*6+1:ii*6));
    end

    % Rotate the path and divide two directions
    [x,y] = centreBox(x,y);    % == make each image in the center of the plot ==
    minX = nanmin(x);
    maxX = nanmax(x);
    subn=x>minX+0.25*(maxX-minX) & x<maxX-0.25*(maxX-minX);
    subx=x(subn);
    suby=y(subn);
    slope=polyfit(subx,suby,1); 
    tAngle=atan(slope(1));
    [posx,posy] = rotatePath(x,y,-1*tAngle);
    post=t;
    % limit=0.15;
    [runData1,runData2] = locateRuns_cz(posx,posy,post,limit);
    
    % add a lap in left direction if necessary
    if runData1(1).dir==1 % the first lap is in right direction
        run0=struct('n',1,'x',posx(1),'y',posy(1),'t',post(1),'dir',2);
        runData1=[run0;runData1];
        runData2=[run0;runData2];
    end
    % ++++++++++++++++ change runData and runData2 +++++++++++++++++++
    runData{aa,1}=runData2;
    nlaps=length(runData{aa,1});
    
    for nl=1:nlaps
        winstart=ceil((runData{aa,1}(nl).n(1)-1)/6)+1;
        winend=floor(runData{aa,1}(nl).n(end)/6);
        runData{aa,1}(nl).nwin=winstart:winend;
    end
end

%% get CSC files
for bb=1:neeg
    S_abs0=[];
    S_rlt0=[];
    S_nwin=[];
    % Get the S_rlt and S_abs for all the sessions
    for aa=1:ndirs  % ndirs is the number of sessions
        filnam = strcat(dataloc{aa},dir_delim,'CSC',num2str(csclist(bb)),'.ncs');
        if(exist(filnam)~=2)
            continue
        end
        [sample, tt, tt_raw, Fs]=loadCSC_new_cz(filnam); % sample unit:uv, tt and tt_raw unit: s
        % gamma_subsample{aa,bb}=ffthighpass(sample, Fs, 20, 25);
        
        % get spectrogram matrix of eeg
        param = struct('tapers',[3,5], 'Fs',Fs, 'pad', padval); % setup as in Ahmed and Mehta 2012: 'tapers',[2 3], 'Fs',Fs
        
        [S_abs{aa,bb}, T, F] = mtspecgramc(sample', [.2 .2], param); %200 ms bins with no overlap, as in Ahmed and Mehta 2012
        S_sum=sum(S_abs{aa,bb},2);
        for nn=1:size(S_abs{aa,bb},1)
            S_rlt{aa,bb}(nn,:)=S_abs{aa,bb}(nn,:)./S_sum(nn);
        end
        S_abs0=[S_abs0;S_abs{aa,bb}];
        S_rlt0=[S_rlt0;S_rlt{aa,bb}];
        S_nwin=[S_nwin;size(S_abs{aa,bb},1)];
    end
    
    % Get the Zscored Sz_rlt and Sz_abs across all the sessions
    Sz_abs0=zscore(S_abs0);
    Sz_rlt0=zscore(S_rlt0);
    n_temp=0;
    for aa=1:ndirs  % ndirs is the number of sessions
        Sz_abs{aa,bb}=Sz_abs0(n_temp+1:n_temp+S_nwin(aa),:);
        Sz_rlt{aa,bb}=Sz_rlt0(n_temp+1:n_temp+S_nwin(aa),:);
        n_temp=n_temp+S_nwin(aa);
    end
    
    % ++++++++++++++++ change Sz_abs and Sz_rlt +++++++++++++++++++
    Sz=Sz_rlt;
    
    nf=length(F);
    minf=find(F<plot_fs(1));
    minf=minf(end);
    maxf=find(F>plot_fs(2));
    maxf=maxf(1);
    nf2=maxf-minf+1;
    
    bincount_left{bb}=zeros(1,nbin);
    spechist_left{bb}=zeros(nf, nbin);
    spechist_left_z{bb}=zeros(nf2, nbin);
    bincount_right{bb}=zeros(1,nbin);
    spechist_right{bb}=zeros(nf, nbin);
    spechist_right_z{bb}=zeros(nf2, nbin);
    bincount_both{bb}=zeros(1,nbin);
    spechist_both{bb}=zeros(nf, nbin);
    spechist_both_z{bb}=zeros(nf2, nbin);
    
    for aa=1:ndirs
        % find average spectrogram for each velocity bin
        nWin=min(size(S_abs{aa,bb},1), length(subvel{aa}));
        nlaps=length(runData{aa,1});
        for nl=1:nlaps
            if isempty(find(badlap{aa}==nl))   % exclude the badlaps
                if runData{aa,1}(nl).dir==1  % find good right laps
                    for cc=runData{aa,1}(nl).nwin
                        if cc<nWin
                            bin=find(vr>=subvel{aa}(cc), 1, 'first');
                            if(~isempty(bin))
                                bin=max(1,bin-1);
                                spechist_right{bb}(:,bin) = spechist_right{bb}(:,bin)+(Sz{aa,bb}(cc,:)');
                                bincount_right{bb}(bin)=bincount_right{bb}(bin) + 1;
                                spechist_both{bb}(:,bin) = spechist_both{bb}(:,bin)+(Sz{aa,bb}(cc,:)');
                                bincount_both{bb}(bin)=bincount_both{bb}(bin) + 1;
                            end
                        end
                    end
                elseif runData{aa,1}(nl).dir==2 % find good left laps
                    for cc=runData{aa,1}(nl).nwin
                        if cc<nWin
                            bin=find(vr>=subvel{aa}(cc), 1, 'first');
                            if(~isempty(bin))
                                bin=max(1,bin-1);
                                spechist_left{bb}(:,bin) = spechist_left{bb}(:,bin)+(Sz{aa,bb}(cc,:)');
                                bincount_left{bb}(bin)=bincount_left{bb}(bin) + 1;
                                spechist_both{bb}(:,bin) = spechist_both{bb}(:,bin)+(Sz{aa,bb}(cc,:)');
                                bincount_both{bb}(bin)=bincount_both{bb}(bin) + 1;
                            end
                        end
                    end
                end
            end
        end
    end
        
    spechist_left0=zeros(nf, nbin);
    spechist_right0=zeros(nf, nbin);
    spechist_both0=zeros(nf, nbin);
    for bin=1:nbin
        if bincount_left{bb}(bin)>0
            spechist_left0(:,bin)= spechist_left{bb}(:,bin)./bincount_left{bb}(bin);
        else
            spechist_left0(:,bin)= spechist_left{bb}(:,bin);
        end
        if bincount_right{bb}(bin)>0
            spechist_right0(:,bin)= spechist_right{bb}(:,bin)./bincount_right{bb}(bin);
        else
            spechist_right0(:,bin)= spechist_right{bb}(:,bin);
        end
        if bincount_both{bb}(bin)>0
            spechist_both0(:,bin)= spechist_both{bb}(:,bin)./bincount_both{bb}(bin);
        else
            spechist_both0(:,bin)= spechist_both{bb}(:,bin);
        end
    end
    spechist_left{bb} = spechist_left0(minf:maxf,:);
    spechist_right{bb} = spechist_right0(minf:maxf,:);
    spechist_both{bb} = spechist_both0(minf:maxf,:);

    % do z-score on the whole frequency band
%     spechist_left_z{bb} = zscore(spechist_left0);
%     spechist_left_z{bb}=spechist_left_z{bb}(minf:maxf,:);
%     spechist_right_z{bb} = zscore(spechist_right0);
%     spechist_right_z{bb}=spechist_right_z{bb}(minf:maxf,:);
%     spechist_both_z{bb} = zscore(spechist_both0);
%     spechist_both_z{bb}=spechist_both_z{bb}(minf:maxf,:);

    % do z-score on the [minf:maxf]
    spechist_left_z{bb} = zscore(spechist_left{bb});
    spechist_right_z{bb} = zscore(spechist_right{bb});
    spechist_both_z{bb} = zscore(spechist_both{bb});
        
    % plot left figures
    ffa_left = figure('Units','normalized','Position',[0 0 0.5625 1]);
    uimagesc(vr(1:end-1)+(vr(2)-vr(1))/2, F(minf:maxf), spechist_left_z{bb}(:, :));
    colorbar;
    axis xy; axis tight; ylim(plot_fs); grid minor;
    tit=['Tetrode ',num2str(csclist(bb)),'-Left'];
    title(tit);
    filenam = strcat('tt',num2str(csclist(bb)),'-left-all sessions');
    saveas(ffa_left, strcat(subdir1,dir_delim,filenam,'.png'),'png');
    close all;
    
    ffa_left = figure('Units','normalized','Position',[0 0 0.5625 1]);
    uimagesc(log2(vr(2:end).*(4/3)), log2(F(minf:maxf).*(4/30)), spechist_left_z{bb}(:, :),plotcolor(1,:));
    colorbar;
    axis xy; axis tight; ylim(log2(plot_fs.*(4/30))); grid minor;
    tickx=0:7;
    set(gca, 'XTick',tickx);
    set(gca, 'XTickLabel',2.^tickx*0.75);
    ticky=2:5;
    set(gca, 'YTick',ticky);
    set(gca, 'YTickLabel',2.^ticky*7.5);
    tit=['Tetrode ',num2str(csclist(bb)),'-Left'];
    title(tit);
    filenam = strcat('tt',num2str(csclist(bb)),'-left-log-all sessions');
    saveas(ffa_left, strcat(subdir1,dir_delim,filenam,'.png'),'png');
    close all;
    
    % plot right figures
    ffa_right = figure('Units','normalized','Position',[0 0 0.5625 1]);
    uimagesc(vr(1:end-1)+(vr(2)-vr(1))/2, F(minf:maxf), spechist_right_z{bb}(:, :));
    colorbar;
    axis xy; axis tight; ylim(plot_fs); grid minor;
    tit=['Tetrode ',num2str(csclist(bb)),'-Right'];
    title(tit);
    filenam = strcat('tt',num2str(csclist(bb)),'-right-all sessions');
    saveas(ffa_right, strcat(subdir1,dir_delim,filenam,'.png'),'png');
    close all;
    
    ffa_right = figure('Units','normalized','Position',[0 0 0.5625 1]);
    uimagesc(log2(vr(2:end).*(4/3)), log2(F(minf:maxf).*(4/30)), spechist_right_z{bb}(:, :),plotcolor(2,:));
    colorbar;
    axis xy; axis tight; ylim(log2(plot_fs.*(4/30))); grid minor;
    tickx=0:7;
    set(gca, 'XTick',tickx);
    set(gca, 'XTickLabel',2.^tickx*0.75);
    ticky=2:5;
    set(gca, 'YTick',ticky);
    set(gca, 'YTickLabel',2.^ticky*7.5);
    tit=['Tetrode ',num2str(csclist(bb)),'-Right'];
    title(tit);
    filenam = strcat('tt',num2str(csclist(bb)),'-right-log-all sessions');
    saveas(ffa_right, strcat(subdir1,dir_delim,filenam,'.png'),'png');
    close all;
    
    % plot bothdir figures
    ffa_both = figure('Units','normalized','Position',[0 0 0.5625 1]);
    uimagesc(vr(1:end-1)+(vr(2)-vr(1))/2, F(minf:maxf), spechist_both_z{bb}(:, :));
    colorbar;
    axis xy; axis tight; ylim(plot_fs); grid minor;
    tit=['Tetrode ',num2str(csclist(bb)),'-Both Dirs'];
    title(tit);
    filenam = strcat('tt',num2str(csclist(bb)),'-both-all sessions');
    saveas(ffa_both, strcat(subdir1,dir_delim,filenam,'.png'),'png');
    close all;
    
    ffa_both = figure('Units','normalized','Position',[0 0 0.5625 1]);
    uimagesc(log2(vr(2:end).*(4/3)), log2(F(minf:maxf).*(4/30)), spechist_both_z{bb}(:, :),plotcolor(3,:));
    colorbar;
    axis xy; axis tight; ylim(log2(plot_fs.*(4/30))); grid minor;
    tickx=0:7;
    set(gca, 'XTick',tickx);
    set(gca, 'XTickLabel',2.^tickx*0.75);
    ticky=2:5;
    set(gca, 'YTick',ticky);
    set(gca, 'YTickLabel',2.^ticky*7.5);
    tit=['Tetrode ',num2str(csclist(bb)),'-Both Dirs'];
    title(tit);
    filenam = strcat('tt',num2str(csclist(bb)),'-both-log-all sessions');
    saveas(ffa_both, strcat(subdir1,dir_delim,filenam,'.png'),'png');
    close all;

end

%% Plot figures
for bb=1:neeg
    % Assemble original figures
    I=[];
    filenam = strcat('tt',num2str(csclist(bb)),'-left-all sessions');
    A = imread(strcat(subdir1,dir_delim,filenam,'.png'));
    I=[I A];
    
    filenam = strcat('tt',num2str(csclist(bb)),'-right-all sessions');
    A = imread(strcat(subdir1,dir_delim,filenam,'.png'));
    I=[I A];
    
    filenam = strcat('tt',num2str(csclist(bb)),'-both-all sessions');
    A = imread(strcat(subdir1,dir_delim,filenam,'.png'));
    I=[I A];
    imwrite(I,strcat(subdir2,dir_delim,'Assemble-tt',num2str(csclist(bb)),'-all.png'));
    
    % Assemble log figures
    I=[];
    filenam = strcat('tt',num2str(csclist(bb)),'-left-log-all sessions');
    A = imread(strcat(subdir1,dir_delim,filenam,'.png'));
    I=[I A];
    
    filenam = strcat('tt',num2str(csclist(bb)),'-right-log-all sessions');
    A = imread(strcat(subdir1,dir_delim,filenam,'.png'));
    I=[I A];
    
    filenam = strcat('tt',num2str(csclist(bb)),'-both-log-all sessions');
    A = imread(strcat(subdir1,dir_delim,filenam,'.png'));
    I=[I A];
    imwrite(I,strcat(subdir2,dir_delim,'Assemble-tt',num2str(csclist(bb)),'-all-log.png'));
end