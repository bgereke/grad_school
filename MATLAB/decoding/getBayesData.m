function [out] = getBayesData

%for monkeying with tfr stuff,since its too large to get all of them

file = 'bayes_3_04_01_2_2_15_Pxxu1000_1';
out = {};

cd begin1
load (file,'scores')
for i = 1:length(scores)    
    scores{i,27} = pwd;    
end
out = [out;scores];
clearvars -except out file
cd ../

cd begin2
load (file,'scores')
for i = 1:length(scores)    
    scores{i,27} = pwd;    
end
out = [out;scores];
clearvars -except out file
cd ../

cd begin3
load (file,'scores')
for i = 1:length(scores)    
    scores{i,27} = pwd;    
end
out = [out;scores];
clearvars -except out 
cd ../

check = 1;
for j = 1:size(out,1)
    if isnan(out{check,22})
        out(check,:) = [];
    else
        check = check +1;
    end   
end