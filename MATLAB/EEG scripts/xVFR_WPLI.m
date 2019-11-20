function [xVFR,mxVFR,zxVFR,mVFRsh] = xVFR_WPLI(ImX,vel,vbins,lags)

invh = 0.5;   
% lags = -lag*30:lag*30;
xVFR = zeros(size(ImX,1),length(vbins),length(lags));

for l = 1:length(lags)
    if lags(l) < 0
        xVFR(:,:,l) = P_estimator(ImX(:,1:end+lags(l)),vel(-lags(l)+1:end),vbins,invh);
    elseif lags(l) > 0        
        xVFR(:,:,l) = P_estimator(ImX(:,lags(l)+1:end),vel(1:end-lags(l)),vbins,invh);
    else
        xVFR(:,:,l) = P_estimator(ImX,vel,vbins,invh);
    end
end

nsh = 500;
VFRsh = zeros(size(xVFR,1),size(xVFR,2),nsh);
shiftmin = 60*29;
idx = randsample(size(ImX,2)-2*shiftmin,nsh,true);
for s = 1:nsh
%    percent_done = s/nsh*100
   ImXshift = [ImX(:,idx(s)+shiftmin+1:end) ImX(:,1:idx(s)+shiftmin)];
   VFRsh(:,:,s) =  P_estimator(ImXshift,vel,vbins,invh);
end

[ci] = prctile(VFRsh,[2.5 50 97.5],3);
mVFRsh = ci(:,:,2);
% nanVFR = VFR>ci(:,:,1) & VFR<ci(:,:,3);
mxVFR = xVFR - repmat(ci(:,:,2),1,1,size(xVFR,3));
% sVFR = mVFR;
% sVFR(nanVFR) = nan;    
zxVFR = (xVFR-repmat(mean(VFRsh,3),1,1,size(xVFR,3)))./repmat(std(VFRsh,0,3),1,1,size(xVFR,3));
% Other Functions

% Estimate the mean P for one position value
function p = P_estimator(P,vel,vbin,invh)
delta = repmat(vel',1,size(vbin,2))-repmat(vbin,length(vel),1);
num = P*gaussian_kernel(delta,invh);
den = abs(P)*gaussian_kernel(delta,invh);
p = (num./den).^2;
                                                    
% Gaussian kernel for the rate calculation
function r = gaussian_kernel(x,invh)
% r =  invh/sqrt(2*pi)*exp(-0.5*x.*x*invh^2);
r =  exp(-0.5*x.*x*invh^2);

% function P = P_estimator(ImX,ImXp,pb,kappa)
% % edge-corrected kernel density estimator
% delta = angle(repmat(exp(1j*ImXp'),1,size(pb,2)).*conj(exp(1j*pb)));
% outsum = ImX*vm_kernel(delta,kappa);
% outsumW = abs(ImX)*vm_kernel(delta,kappa);
% P = (outsum./outsumW).^2;