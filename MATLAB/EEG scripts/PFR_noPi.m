function [PFR] = PFR_noPI(x,ref,freqVec,width,Fs,phasebins)

[theta_phase, thetaBP] = thetaphase(x,Fs);
[noPi] = traces2noPi(x,ref,freqVec,Fs,width);
PFR = zeros(size(noPi,1),size(phasebins,2)-1);
noPi = zscore(noPi,0,2);

for p=1:size(phasebins,2)-1    
    pidx = theta_phase >= phasebins(p) & theta_phase < phasebins(p+1); 
    PFR(:,p) = mean(noPi(:,pidx),2); 
end