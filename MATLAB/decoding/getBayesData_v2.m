function [out] = getBayesData_v2

%for monkeying with tfr stuff,since its too large to get all of them

file = 'bayes_3_04_01_2_2_15_Pxxu1000_1_uV_tfr_w3_onwindows';
out = {};

cd begin1
load (file,'scores')
out = [out;scores];
clearvars -except out file
cd ../

cd begin2
load (file,'scores')
out = [out;scores];
clearvars -except out file
cd ../

cd begin3
load (file,'scores')
out = [out;scores];
clearvars -except out 
cd ../
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%side script for adjusting all decoded time columns so that total P(x|n)
%(for all x) totals to one
%NOTE, do not need to do that in future since it is also writen into
%doall_nopaths, but doesnt matter if done twice anyway...

for i = 1:size(out,1)
    sm = nansum(out{i,12}); %adjust so that each time bin, the sum of P(x|n) (for all x) = 1;
    out{i,21} = out{i,12};
    for o = 1:size(out{i,21},2)
        out{i,21}(:,o) = out{i,21}(:,o)/sm(o);
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
out(:,1:11) = [];
out(:,4:9) = [];
out(:,5:9) = [];
out(:,7) = [];
%get rid of everything except error and tfrs, error is col1, tfr col2

check = 1;
for j = 1:size(out,1)
    if isempty(out{check,5})
        out(check,:) = [];
    else
        out{check,7} = pwd;
        check = check +1;
    end   
end

% 
% check = 1;
% for j = 1:size(out,1)
%     if isnan(out{check,1})
%         out(check,:) = [];
%     else
%         out{check,2} = mean(out{check,2},3);
%         check = check +1;
%     end   
% end