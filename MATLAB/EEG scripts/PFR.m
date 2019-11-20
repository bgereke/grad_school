
function [PFR,cycle] = PFR(S,freqVec,width,Fs,pbins)

[theta_phase,wave_phase,sym,waveBP] = thetaphase(S,Fs);
kappa = 50;
delta = angle(exp(1i*repmat(theta_phase',1,length(pbins))).*conj(exp(1i*repmat(pbins,length(theta_phase),1))));
W = exp(kappa*cos(delta));
D = ones(length(freqVec),length(theta_phase))*exp(kappa*cos(delta));
cycle = S'*W./D;
cycle = cycle(1,:);

if pbins(cycle==min(cycle)) > 0
     S = -S;
    [theta_phase,wave_phase,sym,waveBP] = thetaphase(S,Fs);  
end

[TFR] = traces2TFR_rev(S,freqVec,Fs,width);
TFR = zscore(TFR,0,2);
delta = angle(exp(1i*repmat(theta_phase',1,length(pbins))).*conj(exp(1i*repmat(pbins,length(theta_phase),1))));
W = exp(kappa*cos(delta));
D = ones(length(freqVec),length(theta_phase))*exp(kappa*cos(delta));
cycle = S'*W./D;
cycle = cycle(1,:);
PFR = TFR*W./D;