%run directories_bayesdata; using getBayesData_v2 first

for i = 26229
figure
if length(find(~isnan(allscores{i,4}(1,1:8))))>7
n = 3;
    
subplot(n,1,1);
imagesc([.01:.01:.01*size(allscores{i,4},2)],allscores{i,3},allscores{i,4});
hold on; 
plot([.01:.01:.01*size(allscores{i,2},1)],allscores{i,3}(floor(allscores{i,2})),'r','linewidth',2);
plot([0 .13], [200 200],'k','linewidth',2);
title(i)

% subplot(2,1,2);
% imagesc([.01:.01:.01*size(allscores{i,4},2)],allscores{i,3},allscores{i,1});
% hold on; 
% plot([.01:.01:.01*size(allscores{i,2},1)],allscores{i,3}(floor(allscores{i,2})),'r','linewidth',2);
% plot([0 .13], [200 200],'k','linewidth',2);
% title(i)

% subplot(n,1,2);
 xx = 0:1/2000:(size(allscores{i,5}(:,100:end-100),2)-1)/2000;
% plot(xx,allscores{i,6}(100:end-100,:));
% xlim([min(xx) max(xx)])

subplot(n,1,2);
plot(xx,fftlowpass(mean(allscores{i,6}(100:end-100,:),2),2000,240,250));
xlim([min(xx) max(xx)])

subplot(n,1,3);
yy = 100:-2:25;
imagesc(xx,yy,flipud(allscores{i,5}(:,100:end-100)));
set(gca,'YDir','normal')

%if want to compare the tfr stored to other versions:
%tfr = TFR_frequency_band(mean(allscores{i,6},2),2000,5,25,100);
%figure; imagesc(xx,yy,flipud(tfr);set(gca,'YDir','normal');
  
% pause
% close all
end
end  