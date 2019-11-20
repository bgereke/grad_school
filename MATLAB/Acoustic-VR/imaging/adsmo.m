function [xsm] = adsmo(x,t,lambdax,lambdat)
%adaptive kernel smoothing
dx = repmat(x,1,length(x)) - repmat(x',length(x),1);
smx = exp(-(dx.^2)/lambdax);
smt = zeros(length(x));
for i = 1:length(x) 
    smt(i,:) = exp(-((t-t(i)).^2)/(lambdat));
end
smxt = smx.*smt;
smxt = smxt./repmat(sum(smxt,2),1,length(x));
xsm = smxt*x;


