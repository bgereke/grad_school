function [decays] = expfit(TFRcounts,powbins)

decays = zeros(1,size(TFRcounts,1));

for i = 1:size(TFRcounts,1)
    f = fit(powbins',TFRcounts(i,:)','exp1');
    decays(i) = f.b;
end
    
decays;