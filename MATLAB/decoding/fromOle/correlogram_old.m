function [Cxy,lag] = correlogram(x1,x2,maxlag);

s1(x1) = 1; 
s2(x2) = 1; 

if length(s1) > length(s2)
    s1 = s1(1:length(s2));
else
    s2 = s2(1:length(s1));
end
size(s1)
size(s2)
[Cxy,lag] = xcorr(s1,s2,maxlag,'coeff'); 
