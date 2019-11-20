function Dw = prewhitening(D,Order)
%
% Perform prewhitening using an autoregressive approach where
% filter coefficient are estimated by Yule-Walker approach
%

a = aryule(D,Order);
Dw = filter(a,1,D);