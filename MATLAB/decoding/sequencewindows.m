function [res] = sequencewindows(fields,windows,minspikes)

%minspikes = 3;
res = getthetawindowspikes(fields,windows,minspikes);

for i = 1:size(res,1)   
    res{i,3} = 0; %give default 0 score for all windows
    a = allsequences(res{i,1}); %get all possible seuqnces
    a(:,:,1)
    res{i,4} = a(:,:,1); %default best is the first sequence
    for j = 1:size(a,3)
        score = (scorer(a(:,:,j)));
        if abs(score) > abs(res{i,3})
           res{i,3} = score; 
           res{i,4} = a(:,:,j);
        end
    end
end

%get shuffled score distributions via method 1: shuffle all the COMs to
%existing spike times
numshuffles = 300;
for i = 1:size(res,1)
    
    a = allsequences(res{i,1}); %get all possible sequences for the window
  
%     if size(a,3) > 20
%         pause
%     end
    scores = zeros(1,numshuffles); %reset each shuffled score to 0
    for j = 1:numshuffles
        numspikes = size(a(:,1,1));
        r = randperm(numspikes); %shuffle row numbers
        tempscore = 0;
        for k = 1:size(a,3) 
            
            score = scorer([a(:,1,k) a(r,2,k)]);
            if abs(score) > abs(tempscore)
               tempscore = score;
               scores(j) = score;
            end
        end
    end
    res{i,5} = scores;
    res{i,6} = prctile(scores,5);
    res{i,7} = prctile(scores,95);
end

%get shuffled score distributions via method2: reassign each cell a new COM
for i = 1:size(res,1)
    
    scores = zeros(1,numshuffles);
    for j = 1:numshuffles
        a = allsequencesshuffled(res{i,1});
        %find highest score of the shuffled sequences
        tempscore = 0;
        for k = 1:size(a,3)
            score = scorer([a(:,1,k) a(:,2,k)]);
            if abs(score) > abs(tempscore)
               tempscore = score;
               scores(j) = score;
            end
        end
    end
    res{i,8} = scores;
    res{i,9} = prctile(scores,5);
    res{i,10} = prctile(scores,95);
end

for i = 1:size(res,1)
res{i,11} = ((res{i,3} <= res{i,6} & res{i,3} <= res{i,9})| (res{i,3} >= res{i,7} & res{i,3} >= res{i,10})) & res{i,3} ~= 0;
end