function [VFR,mVFR,zVFR,sVFR] = VFR_WPLI(ImX,vel,vbins)

invh = 0.3;    
VFR = P_estimator(ImX,vel,vbins,invh);
VFR = atanh(VFR);
nsh = 500;
VFRsh = zeros(size(VFR,1),size(VFR,2),nsh);
shiftmin = 60*29;
idx = randsample(size(ImX,2)-2*shiftmin,nsh,true);
for s = 1:nsh
%    percent_done = s/nsh*100
   ImXshift = [ImX(:,idx(s)+shiftmin+1:end) ImX(:,1:idx(s)+shiftmin)];
   VFRsh(:,:,s) =  P_estimator(ImXshift,vel,vbins,invh);
end

VFRsh = atanh(VFRsh);
[ci] = prctile(VFRsh,[2.5 50 97.5],3);
nanVFR = VFR>ci(:,:,1) & VFR<ci(:,:,3);
mVFR = (VFR - ci(:,:,2));
sVFR = mVFR;
sVFR(nanVFR) = nan;    
zVFR = (VFR-mean(VFRsh,3))./std(VFRsh,0,3);
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