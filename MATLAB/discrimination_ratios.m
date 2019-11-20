%LFPmaps(isnan(LFPmaps)) = 0;
% LFPmaps = (LFPmaps-min(min(LFPmaps(LFPmaps ~= 0))))/(max(max(LFPmaps))-min(min(LFPmaps(LFPmaps ~= 0))));
% LFPmaps(LFPmaps < 0) = nan;
LFPmaps = reshape(tiedrank(reshape(LFPmaps,1,60*60)),60,60);
novel = reshape(LFPmaps(:,31:60),1,60*30);
novel(isnan(novel)) = [];
familiar = reshape(LFPmaps(:,1:30),1,60*30);
familiar(isnan(familiar)) = [];
dis_idx = median(novel)/(median(novel)+median(familiar))


novel = reshape(LFPmaps(1:30,1:30),1,30*30);
novel(isnan(novel)) = [];
familiar = reshape(LFPmaps(31:60,31:60),1,30*30);
familiar(isnan(familiar)) = [];
dis_idx = median(novel)/(median(novel)+median(familiar))