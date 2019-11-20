function [PFR] = PFR_ImX(x,y,freqVec,width,Fs,phasebins)

[theta_phase, thetaBP] = thetaphase(x,Fs);
[ImX] = traces2ImX(x,y,freqVec,Fs,width);
PFR = zeros(size(ImX,1),size(phasebins,2)-1);
ImX = abs(ImX);
ImX = zscore(ImX,0,2);

for p=1:size(phasebins,2)-1    
    pidx = theta_phase >= phasebins(p) & theta_phase < phasebins(p+1); 
    PFR(:,p) = mean(ImX(:,pidx),2); 
end