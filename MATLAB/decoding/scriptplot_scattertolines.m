%for a bin
start = 0;
stop = 720;
size = 5;
shift = 5;
numbins = floor((stop - start) / shift);

xx = [Tphase;Tphase+360]; %put data here...
yy = [e;e]; 
% xx(ss>100) = [];
% ss(ss>100) = [];

means = [];    
xxx = [];
errs = [];

startbin = start;
stopbin = start+size;
for i = 1:numbins
    ind = find(xx >= startbin & xx < stopbin);
    means(i) = nanmean(yy(ind));
    errs(i) = 1.96*std(yy(ind))/sqrt(length(ind));
    xxx(i) = startbin;
    startbin = startbin+shift;
    stopbin = stopbin+shift;
    
%     for boxplottin
%     ind = find(xx >= startbin & xx < stopbin);
%     meany = [meany;(yy(ind))];
%         for j = 1:length(ind)
%             xxx = [xxx; startbin];
%         end
%     startbin = startbin+size;
%     stopbin = startbin+size;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for plotting power by error size

% start = -200;
% stop = 210;
% size = 5;
% numbins = floor((stop - start) / size);
% 
% meanf= [];  
% means = []; 
% xxx = [];
% errf = [];
% errs = [];
% 
% startbin = start;
% stopbin = start+size;
% for i = 1:numbins
%     ind = find(yy >= startbin & yy < stopbin);
%     meanf(i) = nanmean(ff(ind));
%     means(i) = nanmean(ss(ind));
%     errs(i) = std(ss(ind))/sqrt(length(ind));
%     errf(i) = std(ff(ind))/sqrt(length(ind));
%     xxx(i) = startbin;
%     startbin = startbin+size;
%     stopbin = startbin+size;
% end
% figure; 
% subplot(2,1,1); plot(xxx,meanf,'k','linewidth',2);hold on;plot(xxx,meanf+errf,'k');plot(xxx,meanf-errf,'k');
% subplot(2,1,2); plot(xxx,means,'k','linewidth',2);hold on;plot(xxx,means+errs,'k');plot(xxx,means-errs,'k');