function [VFR] = VFR_ImX(x,y,freqVec,width,Fs,vel,vbins)

invh = 1;
[ImX] = traces2ImX(x,y,freqVec,Fs,width);
ImX = zscore(abs(ImX),0,2);
VFR = zeros(size(ImX,1),size(vbins,2));

for v=1:size(vbins,2)  
    
    VFR(:,v) = P_estimator(ImX,vel,vbins(v),invh);
    
end

% Other Functions

% Estimate the mean P for one position value
function p = P_estimator(P,vel,vbin,invh)
% edge-corrected kernel density estimator
num = P*gaussian_kernel((vel'-vbin),invh);
den = sum(gaussian_kernel((vel'-vbin),invh));
p = num/den;
                                                    
% Gaussian kernel for the rate calculation
function r = gaussian_kernel(x,invh)
r =  invh/sqrt(2*pi)*exp(-0.5*x.*x*invh^2);