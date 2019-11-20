function [data] = load_gimt_data(mice)

img_text = 'on';

fid = fopen(mice,'r');
if fid == -1
    msgbox('Could not open the input file! Make sure the filname and path are correct.','ERROR');
end

% Get the mice to be analyzed from input file
fid = fopen(mice,'r');
ii = 1;
while ~feof(fid)
    str = fgetl(fid);
    if ~strcmp(str(end),'\')
        str = strcat(str,'\');
    end
    Mice(ii) = {str};
    ii = ii+1;
end
numMice = ii-1;

T = 0:10:700;
Vel = 0:0.5:50;
TGIM = cell(3,1);
VGIM = [];
Mouse = [];
TV = cell(3,1);

for ii = 1:numMice
    disp(sprintf('%s%s','Reading data for: ',Mice{ii}));
    
    % find and open the dates.txt file
    dates_txt = strcat(Mice{ii},'dates.txt');
    fid = fopen(dates_txt,'r');
    jj = 1;
    while ~feof(fid)
        str = fgetl(fid);
        if ~strcmp(str(end),'\')
            str = strcat(str,'\');
        end
        Dates(jj) = {str};
        jj = jj+1;
    end
    numDates = jj-1;
    
    % Perform all analyses for each date
    for jj=1:numDates
        disp(Dates{jj});
        newpath = strcat(Mice{ii},Dates{jj});
        cd(newpath);
        
        load(strcat(cd,'\GIMt\data.mat'));
        
        for kk = 1:3
            TGIM{kk} = cat(3,TGIM{kk},[zGIM{kk} nan(size(GIM{kk},1),length(T)-size(GIM{kk},2))]);
            TV{kk} = cat(2,TV{kk},[V{kk} nan(1,length(T)-size(GIM{kk},2))]);
        end
        
        load(strcat(cd,'\VFR_GIM\data.mat'));   
        VGIM = cat(3,VGIM,[zGIM nan(size(zGIM,1),length(Vel)-size(zGIM,2))]);

        Mouse = [Mouse;ii];        
    end    
end

close all

zGIM1 = nanmean(TGIM{1}(:,T<=600,:),3);zGIM2 = nanmean(TGIM{2}(:,T<=600,:),3);zGIM3 = nanmean(TGIM{3}(:,T<=600,:),3);
sp = 60;
c = max(max(abs([zGIM1(freqVec>20,:) zGIM2(freqVec>20,:) zGIM3(freqVec>20,:)])));
subplot(2,1,1)
imagesc(T(T<=600),freqVec,zGIM1,[-c c]);hold on
tm = 600;
plot([tm+0.5*sp tm+0.5*sp],[0 max(freqVec)+1],'-k')
imagesc(T(T<=600)+tm+sp,freqVec,zGIM2,[-c c]);hold on
tm = max(T(T<=600)+tm+sp);
plot([tm+0.5*sp tm+0.5*sp],[0 max(freqVec)+1],'-k')
imagesc(T(T<=600)+tm+sp,freqVec,zGIM3,[-c c]);
axis xy;hold off
xlabel('time (sec)');ylabel('frequency (Hz)');
xlim([0 max(T(T<=600)+tm+sp)])
cbfreeze
freezeColors
cblabel('zGIM')

numperms = 2000;
pzGIM1 = zeros(size(zGIM1,1),size(zGIM1,2),numperms);pzGIM2 = pzGIM1;pzGIM3 = pzGIM1;
TGIM{1} = TGIM{1}(:,T<=600,:);
TGIM{2} = TGIM{2}(:,T<=600,:);
TGIM{3} = TGIM{3}(:,T<=600,:);
for p = 1:numperms
    start1 = randsample(1:size(TGIM{1},3),size(TGIM{1},3),'true');
    start2 = randsample(1:size(TGIM{1},3),size(TGIM{2},3),'true');
    start3 = randsample(1:size(TGIM{1},3),size(TGIM{3},3),'true');

    pzGIM1(:,:,p) = nanmean(TGIM{1}(:,:,start1),3);
    pzGIM2(:,:,p) = nanmean(TGIM{2}(:,:,start2),3);
    pzGIM3(:,:,p) = nanmean(TGIM{3}(:,:,start3),3);
end

zzGIM1 = (zGIM1)./nanstd(pzGIM1,0,3);
zzGIM2 = (zGIM2)./nanstd(pzGIM2,0,3);
zzGIM3 = (zGIM3)./nanstd(pzGIM3,0,3);

sp = 60;th = 3;
c = max(max(abs([zzGIM1(freqVec>20,:) zzGIM2(freqVec>20,:) zzGIM3(freqVec>20,:)])));
subplot(2,1,2)
imagesc(T(T<=600),freqVec,zzGIM1,[-c c]);hold on
contour(T(T<=600),freqVec,zzGIM1,[-th th],'-k','LineWidth',2)
tm = 600;
plot([tm+0.5*sp tm+0.5*sp],[0 max(freqVec)+1],'-k')
imagesc(T(T<=600)+tm+sp,freqVec,zzGIM2,[-c c]);hold on
contour(T(T<=600)+tm+sp,freqVec,zzGIM2,[-th th],'-k','LineWidth',2)
tm = max(T(T<=600)+tm+sp);
plot([tm+0.5*sp tm+0.5*sp],[0 max(freqVec)+1],'-k')
imagesc(T(T<=600)+tm+sp,freqVec,zzGIM3,[-c c]);
contour(T(T<=600)+tm+sp,freqVec,zzGIM3,[-th th],'-k','LineWidth',2)
axis xy;hold off
xlabel('time (sec)');ylabel('frequency (Hz)');
xlim([0 max(T(T<=600)+tm+sp)])
cbfreeze
freezeColors
cblabel('zzGIM')

figure
vmax = 25;
zGIM1 = nanmean(VGIM(:,Vel<=vmax,:),3);
c = max(max(abs(zGIM1)));
subplot(2,1,1)
imagesc(Vel(Vel<=vmax),freqVec,zGIM1,[-c c]);hold on
axis xy;axis square;hold off;%colormap hot
xlabel('runnning speed (cm/sec)');ylabel('frequency (Hz)');
cbfreeze
freezeColors
cblabel('zGIM')

numperms = 2000;
pzGIM1 = zeros(size(zGIM1,1),size(zGIM1,2),numperms);
VGIM = VGIM(:,Vel<=vmax,:);
for p = 1:numperms
    start1 = randsample(1:size(VGIM,3),size(VGIM,3),'true');

    pzGIM1(:,:,p) = nanmean(VGIM(:,:,start1),3);
end

zzGIM1 = (zGIM1)./nanstd(pzGIM1,0,3);

th = 3;
c = max(max(abs(zzGIM1(freqVec>20,:))));
subplot(2,1,2)
imagesc(Vel(Vel<=vmax),freqVec,zzGIM1,[-c c]);hold on
contour(Vel(Vel<=vmax),freqVec,zzGIM1,[-th th],'-k','LineWidth',2)
axis xy;axis square;hold off;colormap default
xlabel('runnning speed (cm/sec)');ylabel('frequency (Hz)');
cbfreeze
freezeColors
cblabel('zzGIM')


% figImage = sprintf('%s%s%s%s',strcat('GIMt','\GIMt'),'.fig');
% bmpImage = sprintf('%s%s%s%s',strcat('GIMt','\GIMt'),'.bmp');
% saveas(gcf,figImage,'fig');
% saveas(gcf,bmpImage,'bmp');
% 
% data.zGIM = zGIM; data.freqVec = freqVec;
% filename = sprintf('%s%s%s%s',strcat('GIMt','\'),'data.mat');
% save(filename,'-struct','data');


