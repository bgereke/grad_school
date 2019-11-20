function [fit1,m1] = localfit(t,x,tspan)

tspan = tspan/2;
t = reshape(t,length(t),1);  
fit1 = zeros(size(x));
m1 = nan(size(x));
x = reshape(x,length(x),1);
for i = 1:length(t)
    support = t(abs(t(i)-t)<=tspan);
    targets = x(abs(t(i)-t)<=tspan);
    if isnan(x(i))
        fit1(i) = nan;
    elseif sum(~isnan(targets))==1
        fit1(i) = x(i);
    else
        tar = targets(~isnan(targets));
        sup = support(~isnan(targets));
        %solve regression using QR decomposition
        [Q,R] = qr([ones(size(sup)),sup],0);
        b = R\(Q'*tar);        
        fit1(i) = b(2)*t(i)+b(1);
        m1(i) = b(2);
    end
end
