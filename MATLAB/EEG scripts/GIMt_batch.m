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

function [GIMt] = GIMt_batch(inFile,sigma,numperms)

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

dirInfo = dir(cd);
found = 0;
for kk=1:size(dirInfo,1)
    if dirInfo(kk).isdir
        if strcmp(dirInfo(kk).name,strcat('GIMt','\'))
            found = 1;
        end
    end
end
if found==0
    mkdir(cd,strcat('GIMt','\'));
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

zGIM = cell(numsessions,1);
data.GIM = cell(numsessions,1);
data.t = cell(numsessions,1);
data.V = cell(numsessions,1);

for ii = 1:numsessions
    
    disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));    
    
%     file = strcat(sessions{ii},'vt1.nvt');
%     [t_t, x_t, y_t] = Nlx2MatVT(file,fieldSelection,extractHeader,extractMode);
%     ind = find(x_t == 0);
%     t_t(ind) = [];
%     x_t(ind) = [];
%     y_t(ind) = [];
%     % Do median filtering to supress off-values
%     %[x, y, t] = medianFilter(x, y, t);
%     % Smoothing the position samples with a simple mean filter
%     for cc=8:length(x_t)-7
%         x_t(cc) = nanmean(x_t(cc-7:cc+7));
%         y_t(cc) = nanmean(y_t(cc-7:cc+7));
%     end
%     y_t = -y_t + max(y_t)-min(y_t) + 2*max(y_t); % reflects posisitons so they are consistent with orientation in recording room
%     [hx,xc] = hist(x_t,50); [hy,yc] = hist(y_t,50);
%     th = 50;
%     xmin = xc(find(hx>th,1,'first'));
%     xmax = xc(find(hx>th,1,'last'));
%     ymin = yc(find(hy>th,1,'first'));
%     ymax = yc(find(hy>th,1,'last'));
%     
%     % Adjust input data
%     scalex = 100/(xmax-xmin);
%     scaley = 100/(ymax-ymin);
%     x_t = x_t * scalex;
%     y_t = y_t * scaley;
%     xmin = xmin * scalex;
%     ymin = ymin * scaley;
%     xmax = xmax * scalex;
%     ymax = ymax * scaley;
%     xcen = xmin+(xmax-xmin)/2;
%     ycen = ymin+(ymax-ymin)/2;
%     x_t = x_t - xcen;
%     y_t = y_t - ycen;
%     
%     % Convert timestamps to seconds
%     t_t = t_t/1000000;
%     
%     vel = zeros(length(x_t),1);
%     for v = 2:1:(length(x_t)-1)
%         vel(v) = sqrt((x_t(v+1)-x_t(v-1))^2 + (y_t(v+1)-y_t(v-1))^2)/(t_t(v+1)-t_t(v-1));
%     end
%     vel(end) = sqrt((x_t(end)-x_t(end-1))^2 + (y_t(end)-y_t(end-1))^2)/(t_t(length(x_t))-t_t(length(x_t)-1));
%     vel(1) = vel(2);
%     vel(vel>=40) = 0.5*(vel(circshift((vel>=40),[-3,0]))+vel(circshift((vel>=40),[3,0])));
%     vel = smooth(vel,15);
%     
%     t_t  = t_t - min(t_t); 
    
    load(strcat(sessions{ii},'CS\CS.mat'));
%     CS(:,ind,:,:) = [];
    t_t = time-min(time);
    t = 0:0.5*sigma:max(t_t);
    
    start = linspace(sigma+0.5*sigma,max(t_t)-sigma-0.5*sigma,numperms);
    for p = 1:numperms
        [~,start(p)] = min((t_t - start(p)).^2);
    end
  
    delta = repmat(t_t,1,length(t))-repmat(t,length(t_t),1);
    W = exp(-0.5*delta.*delta/sigma^2);
    D = ones(length(freqVec),length(t_t))*exp(-0.5*delta.*delta/sigma^2);
    MC = zeros(numchannels,numchannels,length(freqVec),length(t));
    pMC = zeros(numchannels,numchannels,length(freqVec),length(t),numperms);
    for jj = 1:numchannels
        for kk = jj:numchannels      
            MC(jj,kk,:,:) = CS(:,:,jj,kk)*W./D;
            for p = 1:numperms
                pMC(jj,kk,:,:,p) = [CS(:,start(p):end,jj,kk) CS(:,1:start(p)-1,jj,kk)]*W./D;
            end
        end
    end
    clear W D
    
    GIM = zeros(length(freqVec),length(t));
    pGIM = zeros(length(freqVec),length(t),numperms);
    
    for f = 1:length(freqVec)
        for it  = 1:length(t)
            mc = MC(:,:,f,it) - tril(MC(:,:,f,it),-1) + tril(MC(:,:,f,it)',-1);
            mc = mc./mean(mean(abs(mc)));
            GIM(f,it) = 0.5*trace(inv(real(mc))*imag(mc)*inv(real(mc))*imag(mc)');
        end
        for p = 1:numperms
            for it  = 1:length(t)
                mc = pMC(:,:,f,it,p) - tril(pMC(:,:,f,it,p),-1) + tril(pMC(:,:,f,it,p)',-1);
                mc = mc./mean(mean(abs(mc)));
                pGIM(f,it,p) = 0.5*trace(inv(real(mc))*imag(mc)*inv(real(mc))*imag(mc)');
            end
        end
    end

    data.V{ii} = vel'*exp(-0.5*delta.*delta/sigma^2)./(ones(1,length(t_t))*exp(-0.5*delta.*delta/sigma^2));
   
    zGIM{ii} = (GIM - mean(pGIM,3))./std(pGIM,0,3); 
    data.GIM{ii} = GIM;
    data.t{ii} = t;
end

close all

sp = 60;
maxcol = [max(max(data.GIM{1}(freqVec>6,:))) max(max(data.GIM{2}(freqVec>6,:))) max(max(data.GIM{3}(freqVec>6,:)))];
maxcol = max(maxcol);
mincol = [min(min(data.GIM{1}(freqVec>6,:))) min(min(data.GIM{2}(freqVec>6,:))) min(min(data.GIM{3}(freqVec>6,:)))];
mincol = min(mincol);
subplot(3,1,1)
imagesc(data.t{1},freqVec,data.GIM{1},[mincol maxcol]);hold on
tm = max(data.t{1});
plot([tm+0.5*sp tm+0.5*sp],[0 max(freqVec)+1],'-k')
imagesc(data.t{2}+max(data.t{1})+sp,freqVec,data.GIM{2},[mincol maxcol]);hold on
tm = max(data.t{2}+max(data.t{1})+sp);
plot([tm+0.5*sp tm+0.5*sp],[0 max(freqVec)+1],'-k')
imagesc(data.t{3}+max(data.t{1})+max(data.t{2})+2*sp,freqVec,data.GIM{3},[mincol maxcol]);
axis xy;colormap hot;hold off
xlabel('time (sec)');ylabel('frequency (Hz)');
xlim([0 max(data.t{3}+max(data.t{1})+max(data.t{2})+2*sp)])
cbfreeze
freezeColors
cblabel('GIM')


col = [max(max(abs(zGIM{1}(freqVec>6,:)))) max(max(abs(zGIM{2}(freqVec>6,:)))) max(max(abs(zGIM{3}(freqVec>6,:))))];
col = max(col);
subplot(3,1,2)
imagesc(data.t{1},freqVec,zGIM{1},[-col col]);axis xy;hold on
tm = max(data.t{1});
plot([tm+0.5*sp tm+0.5*sp],[0 max(freqVec)+1],'-k')
imagesc(data.t{2}+max(data.t{1})+sp,freqVec,zGIM{2},[-col col]);axis xy;
tm = max(data.t{2}+max(data.t{1})+sp);
plot([tm+0.5*sp tm+0.5*sp],[0 max(freqVec)+1],'-k')
imagesc(data.t{3}+max(data.t{1})+max(data.t{2})+2*sp,freqVec,zGIM{3},[-col col]);
axis xy;colormap jet;colorbar
ylabel('frequency (Hz)');
xlim([0 max(data.t{3}+max(data.t{1})+max(data.t{2})+2*sp)])
cblabel('z-score from null')

subplot(3,1,3)
plot(data.t{1},data.V{1},'-k','LineWidth',2);hold on
tm = max(data.t{1});
plot([tm+0.5*sp tm+0.5*sp],[0 max([data.V{1} data.V{2} data.V{3}])+1],'-k')
plot(data.t{2}+max(data.t{1})+sp,data.V{2},'-k','LineWidth',2);
tm = max(data.t{2}+max(data.t{1})+sp);
plot([tm+0.5*sp tm+0.5*sp],[0 max([data.V{1} data.V{2} data.V{3}])+1],'-k')
plot(data.t{3}+max(data.t{1})+max(data.t{2})+2*sp,data.V{3},'-k','LineWidth',2);
ylabel('running speed (cm/sec)');
xlim([0 max(data.t{3}+max(data.t{1})+max(data.t{2})+2*sp)])
ylim([0 max([data.V{1} data.V{2} data.V{3}])+1]);
colorbar; %just for plot size

figImage = sprintf('%s%s%s%s',strcat('GIMt','\GIMt'),'.fig');
bmpImage = sprintf('%s%s%s%s',strcat('GIMt','\GIMt'),'.bmp');
saveas(gcf,figImage,'fig');
saveas(gcf,bmpImage,'bmp');

data.zGIM = zGIM; data.freqVec = freqVec;
filename = sprintf('%s%s%s%s',strcat('GIMt','\'),'data.mat');
save(filename,'-struct','data');


