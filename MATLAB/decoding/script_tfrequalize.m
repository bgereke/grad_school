%to average all the tfrs from bayes theta windows since they are not all
%the same size

%note: make sure to run getBayesData_v2 BEFORE running this

%average all the traces from each window
for i = 1:size(allscores,1)
   allscores{i,3} = (mean(allscores{i,3},2))'; 
end


maxx = 0;
for i = 1:size(allscores,1)
    temp = size(allscores{i,2},2);
    if temp > maxx
        maxx = temp;
    end    
end

%transform all so they are the same size, by adding on extra nan
%padding,with all centered by their middle index

for i = 1:size(allscores,1)
    diffsize = maxx - size(allscores{i,2},2);
    %must add half that amount on either side
    if mod(diffsize,2) %returns 1 if odd
        leftadd = floor(diffsize / 2);
        rightadd = ceil(diffsize / 2);
    else
        leftadd = diffsize /2;
        rightadd = diffsize / 2;
    end
    numrows = size(allscores{i,2},1);
    leftadd = nan(numrows,leftadd);
    rightadd = nan(numrows,rightadd);
    allscores{i,2} = [leftadd,allscores{i,2},rightadd]; 
    %now normalize all power values to max = 1
%     maxpow = nanmax(nanmax(allscores{i,2}));
%     allscores{i,2} = allscores{i,2}/maxpow;
    
     %do same for traces;
     
     leftadd (2:end,:) = [];
     rightadd(2:end,:) = [];
     allscores{i,3} = [leftadd,allscores{i,3},rightadd];
end

%now go through and collect the tfrs in one big matrix depending on
%whatever condition


tfrneg = zeros(size(allscores{1,2}));
tfrpos = zeros(size(allscores{1,2}));
tfrzer = zeros(size(allscores{1,2}));
traceneg = zeros(size(allscores{1,3}));
tracepos = zeros(size(allscores{1,3}));
tracezer = zeros(size(allscores{1,3}));

tfrnegcount=zeros(size(allscores{1,2}));
tfrposcount=zeros(size(allscores{1,2}));
tfrzercount=zeros(size(allscores{1,2}));
tracenegcount=zeros(size(allscores{1,3}));
traceposcount=zeros(size(allscores{1,3}));
tracezercount=zeros(size(allscores{1,3}));


for i = 1:size(allscores,1)
    if allscores{i,1} < -3
        tfrnegcount = tfrnegcount + ~isnan(allscores{i,2});
        tracenegcount = tracenegcount + ~isnan(allscores{i,3});
        tfrneg = nansum(cat(3,tfrneg, allscores{i,2}),3);
        allscores{i,2} = nan;
        traceneg = nansum([traceneg;allscores{i,3}],1);
        allscores{i,3} = nan;
    elseif allscores{i,1}>3 
        tfrposcount = tfrposcount + ~isnan(allscores{i,2});
        traceposcount = traceposcount + ~isnan(allscores{i,3});
        tfrpos = nansum(cat(3,tfrpos, allscores{i,2}),3);
        allscores{i,2} = nan;
        tracepos = nansum([tracepos;allscores{i,3}],1);
        allscores{i,3} = nan;
    else
        tfrzercount = tfrzercount + ~isnan(allscores{i,2});
        tracezercount = tracezercount + ~isnan(allscores{i,3});
        tfrzer = nansum(cat(3,tfrzer, allscores{i,2}),3);
        allscores{i,2} = nan;
        tracezer = nansum([tracezer;allscores{i,3}],1);
        allscores{i,3} = nan;
    end
end

avgnegtfr = flipud(tfrneg./tfrnegcount);
avgpostfr = flipud(tfrpos./tfrposcount);
avgzertfr = flipud(tfrzer./tfrzercount);

avgnegtrace = (traceneg./tracenegcount);
avgpostrace = (tracepos./traceposcount);
avgzertrace = (tracezer./tracezercount);

allcount = (tfrnegcount+tfrposcount+tfrzercount);
avgalltfr = avgnegtfr.*tfrnegcount+avgpostfr.*tfrposcount+avgzertfr.*tfrzercount;
avgalltfr = avgalltfr./allcount;


