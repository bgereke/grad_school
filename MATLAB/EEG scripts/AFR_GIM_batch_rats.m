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

function [GIMt] = AFR_GIM_batch_rats(inFile,sigma,numperms)

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
        if strcmp(dirInfo(kk).name,strcat('AFR_MC','\'))
            found = 1;
        end
    end
end
if found==0
    mkdir(cd,strcat('AFR_MC','\'));
end

% % Get position data
% fieldSelection(1) = 1; % Timestamps
% fieldSelection(2) = 1; % Extracted X
% fieldSelection(3) = 1; % Extracted Y
% fieldSelection(4) = 0; % Extracted Angel
% fieldSelection(5) = 0; % Targets
% fieldSelection(6) = 0; % Points
% % Do we return header 1 = Yes, 0 = No.
% extractHeader = 0;
% % 5 different extraction modes, see help file for Nlx2MatVt
% extractMode = 1; % Extract all data

%compute cross-spectra and running speed for each session
fCS = [];
V = [];

for ii = 1:numsessions
    
    disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));
    
%     %compute running speed
%     file = strcat(sessions{ii},'vt1.nvt');
%     [t_t, x_t, y_t] = Nlx2MatVT(file,fieldSelection,extractHeader,extractMode);
%     
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

    [data] = load(strcat(sessions{ii},'CS\CS.mat'));
%     CS(:,ind,:,:) = [];
    freqVec = data.freqVec;
    PCS = imag(data.PCS).*repmat(sign(mean(imag(data.PCS),2)),1,size(data.PCS,2));
%     PCS = abs(TFR);
    PCS = zscore(PCS,0,2);
    fCS = [fCS PCS];clear CS
    V = [V;data.acc];
end
% keyboard
%compute running speed triggered cross spectra
vmin = quantile(V,0);
vmax = quantile(V,1);
vbins = vmin:2:vmax;
pMC = zeros(length(freqVec),length(vbins),numperms);
start = round(linspace(100,length(V)-100,numperms));
delta = repmat(V,1,length(vbins))-repmat(vbins,length(V),1);
W = exp(-0.5*delta.*delta/sigma^2);
D = ones(length(freqVec),length(V))*exp(-0.5*delta.*delta/sigma^2);
% fCS = imag(fCS).*repmat(sign(mean(imag(fCS),2)),1,size(fCS,2));
% fCS = zscore(fCS,0,2);

MC = fCS*W./D;
for p = 1:numperms
    pMC(:,:,p) = [fCS(:,start(p):end) fCS(:,1:start(p)-1)]*W./D;
end

clear W D

% %compute running speed triggered GIM
% GIM = zeros(length(freqVec),length(vbins));
% pGIM = zeros(length(freqVec),length(vbins),numperms);
% MC = permute(MC,[3,4,1,2]); %reorder dimensions to improve loop speed
% pMC = permute(pMC,[3,4,1,2,5]);
% 
% for f = 1:length(freqVec)
%     for iv  = 1:length(vbins)
%         mc = MC(:,:,f,iv) - tril(MC(:,:,f,iv),-1) + tril(MC(:,:,f,iv)',-1);
%         mc = mc./mean(mean(abs(mc)));
%         GIM(f,iv) = 0.5*trace(inv(real(mc))*imag(mc)*inv(real(mc))*imag(mc)');
%     end
%     for p = 1:numperms
%         for iv  = 1:length(vbins)
%             mc = pMC(:,:,f,iv,p) - tril(pMC(:,:,f,iv,p),-1) + tril(pMC(:,:,f,iv,p)',-1);
%             mc = mc./mean(mean(abs(mc)));
%             pGIM(f,iv,p) = 0.5*trace(inv(real(mc))*imag(mc)*inv(real(mc))*imag(mc)');
%         end
%     end
% end

zMC = (MC - mean(pMC,3))./std(pMC,0,3);
data.MC = MC;clear MC

close all

col = max(max(abs(data.MC(freqVec>6,:))));
subplot(2,1,1)
imagesc(vbins,freqVec,data.MC,[-col col]);
axis xy;axis square;colormap jet;hold off
xlabel('acceleration (cm/sec^2)');ylabel('frequency (Hz)');
cbfreeze
freezeColors
cblabel('fCSz')

col = max(max(abs(zMC(freqVec>6,:))));
subplot(2,1,2)
imagesc(vbins,freqVec,zMC,[-col col]);
axis xy;axis square;colormap jet;colorbar
xlabel('acceleration (cm/sec^2)');ylabel('frequency (Hz)');
cblabel('z-score from null')

figImage = sprintf('%s%s%s%s',strcat('AFR_MC','\AFR_MC'),'.fig');
bmpImage = sprintf('%s%s%s%s',strcat('AFR_MC','\AFR_MC'),'.bmp');
saveas(gcf,figImage,'fig');
saveas(gcf,bmpImage,'bmp');

data.zMC = zMC; data.freqVec = freqVec; data.vbins = vbins;
filename = sprintf('%s%s%s%s',strcat('AFR_MC','\'),'data.mat');
save(filename,'-struct','data');
