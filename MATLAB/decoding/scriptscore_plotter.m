% %for plotting:
%figure;
xx = [];
yy =[];
ffz = []; ssz = [];
ff = []; ss = [];
dir = {};
w =[];
for i = 1:size(allscores,1)
    if ~isempty(allscores{i,22}) %&& allscores{i,13}(1,1)<135 && allscores{i,13}(1,1)>0 && allscores{i,11} == 1  
    %hold on;
    %
     x = allscores{i,14}(round(median(allscores{i,13})));
     if x > max(allscores{i,14})/2
         x = x-max(allscores{i,14})/2;
     end
    %
%     if allscores{i,14}(ceil(allscores{i,13}(1,1)),1)> 200
%      x = allscores{i,14}(ceil(allscores{i,13}(1,1)),1)- 200;
%     else 
%      x = allscores{i,14}(ceil(allscores{i,13}(1,1)),1);
%     end
    %plot(x,  allscores{i,22},'ko');
    xx = [xx;x];
    yy = [yy;allscores{i,22}];
     ff = [ff;allscores{i,23}];
    ss = [ss;allscores{i,24}];
    ffz = [ffz;allscores{i,25}];
    ssz = [ssz;allscores{i,26}];
    dir = [dir;allscores{i,27}];
    w = [w;allscores{i,2}];
    end
end
clear i x
 ff = ff(~isnan(yy));
 ss = ss(~isnan(yy));
 ffz = ffz(~isnan(yy));
 ssz = ssz(~isnan(yy));
 xx = xx(~isnan(yy));
 yy = yy(~isnan(yy));
 beta = polyfit(xx,yy,1);
 px = 0:1:200;
 py = beta(2)+beta(1)*px;
 ff(yy>200 | yy<-200) = [];
 ss(yy>200 | yy<-200) = [];
 ffz(yy>200 | yy<-200) = [];
 ssz(yy>200 | yy<-200) = [];
 xx(yy>200 | yy<-200) = [];
 yy(yy>200 | yy<-200) = [];
%  plot(xx,yy,'ko');hold on;plot ([0 200],[0 0],'r','linewidth',2);hold on; plot(px,py,'linewidth',2);
 