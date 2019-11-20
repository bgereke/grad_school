% %for getting the score of each one
actual = []; 
best = [];
allcom = [];
score = 0;
for i = 1:size(scores,1)
%i
    for g = 1:size(scores{i,21},2)
        all = 0; 
        for j = 1:size(scores{i,21},1)
             all = nansum([all  j*scores{i,21}(j,g)]);
        end
        com = all/nansum(scores{i,21}(:,g));
        allcom (1,g) = com;
    end
    
%%%%%%%%%%%%%scoremethod1
    c = 0;
    for m=1:length(allcom)
        if ~isnan(allcom(m))           
            score = score + (scores{i,14}(floor(allcom(m)))-scores{i,14}(floor(scores{i,13}(m)))  );
            c = c+1;
        end        
    end
    scores{i,22} = score/c; 
    
%%%%%%%%%%%%%scoremethod2 
%     if isempty(find(~isnan(allcom), 1))
%         scores{i,22} = nan;
%     else        
%     ind = find(~isnan(allcom));
%     f = scores{i,14}(floor(allcom(ind))) - scores{i,14}(floor(scores{i,13}(ind)));
%     score = nansum(f);  
%     scores{i,22} = score;
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%    
    clear actual best allcom;
    score = 0;
end

