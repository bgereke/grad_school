function [PFR] = PFR_WRFR(x,y,freqVec,width,Fs,phasebins)

[theta_phase, thetaBP] = thetaphase(x,Fs);
[X] = traces2X(x,y,freqVec,Fs,width);
R = abs(X);
R = zscore(R,0,2);
ImX = imag(X);
ImX = abs(ImX);
ImX = zscore(ImX,0,2);
PFR = zeros(size(X,1),size(phasebins,2)-1);

for p=1:size(phasebins,2)-1    
    pidx = theta_phase >= phasebins(p) & theta_phase < phasebins(p+1); 
    PFR(:,p) = mean(ImX(:,pidx),2);
end