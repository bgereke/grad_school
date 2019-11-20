
function [PFR,cycle] = PFR_rats(S,freqVec,width,Fs,pbins,vel)

minr = 10;
[theta_phase,wave_phase,sym,waveBP] = thetaphase(S,Fs);
theta_phase = theta_phase(vel>minr);
kappa = 50;
delta = angle(exp(1i*repmat(theta_phase',1,length(pbins))).*conj(exp(1i*repmat(pbins,length(theta_phase),1))));
W = exp(kappa*cos(delta));
D = ones(length(freqVec),length(theta_phase))*exp(kappa*cos(delta));
cycle = S(vel>minr)'*W./D;
cycle = cycle(1,:)';

[TFR] = traces2TFR_rev(S,freqVec,Fs,width);
TFR = TFR(:,vel>minr);
TFR = zscore(TFR,0,2);

if sum(sum(TFR(freqVec>65&freqVec<100,theta_phase>-pi/2&theta_phase<pi/2))) < 0
     S = -S;
    [theta_phase,wave_phase,sym,waveBP] = thetaphase(S,Fs);  
    theta_phase = theta_phase(vel>minr);
end

delta = angle(exp(1i*repmat(theta_phase',1,length(pbins))).*conj(exp(1i*repmat(pbins,length(theta_phase),1))));
W = exp(kappa*cos(delta));
D = ones(length(freqVec),length(theta_phase))*exp(kappa*cos(delta));
cycle = S(vel>minr)'*W./D;
cycle = cycle(1,:);
PFR = TFR*W./D;