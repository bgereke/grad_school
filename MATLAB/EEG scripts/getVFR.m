function [VFR] = getVFR(TFR,freqVec,width,Fs,vel,vbins)

invh = 0.25;
VFR = zeros(size(TFR,1),size(vbins,2));

for v=1:size(vbins,2)  
    
    VFR(:,v) = P_estimator(TFR,vel,vbins(v),invh);
    
%     vidx = vel >= vbins(v) & vel < vbins(v+1); 
%     outsum   = nansum(ImX(:,vidx),2);      % compute the sum; 
%     outsumW  = nansum(abs(ImX(:,vidx)),2); % normalization of the WPLI
%     VFR(:,v) = abs(outsum)./outsumW;
end

% Other Functions

% Estimate the mean P for one position value
function p = P_estimator(P,vel,vbin,invh)
% edge-corrected kernel density estimator
num = P*gaussian_kernel((vel-vbin),invh);
den = ones(size(P))*gaussian_kernel((vel-vbin),invh);
p = num./den;
                                                    
% Gaussian kernel for the rate calculation
function r = gaussian_kernel(x,invh)
r =  invh/sqrt(2*pi)*exp(-0.5*x.*x*invh^2);
