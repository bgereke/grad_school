function [PFR,mPFR,zPFR,sPFR] = PFR_WPLI(x,y,freqVec,width,Fs,numbins)

%determine which channel has higher theta power
[Px,~] = doPwelch(x,Fs);
[Py,F] = doPwelch(y,Fs);
thetax = sum(Px(F>6 & F<12));
thetay = sum(Py(F>6 & F<12));
if thetax > thetay
    [theta_phase, ~] = thetaphase(x,Fs);
else
    [theta_phase, ~] = thetaphase(y,Fs);
end

[ImX] = traces2ImX(x,y,freqVec,Fs,width);

phasebins = repmat(linspace(-pi,pi,numbins),size(ImX,2),1);
PFR = P_estimator(ImX,theta_phase,phasebins,500); 

nsh = 100;
PFRsh = zeros(size(PFR,1),size(PFR,2),nsh);
idx = reshape(randsample(size(ImX,2),size(ImX,2)*nsh,true),nsh,size(ImX,2));
for s = 1:nsh
%    percent_done = s/nsh*100
   PFRsh(:,:,s) =  P_estimator(ImX(:,idx(s,:)),theta_phase,phasebins,500);
end

[ci] = prctile(PFRsh,[2.5 0.5 97.5],3);
nanPFR = PFR>ci(:,:,1) & PFR<ci(:,:,3);
mPFR = (PFR - ci(:,:,2));
sPFR = mPFR;
sPFR(nanPFR) = nan;
zPFR = (PFR-mean(PFRsh,3))./std(PFRsh,0,3);

function P = P_estimator(ImX,ImXp,pb,kappa)
% edge-corrected kernel density estimator
delta = angle(repmat(exp(1j*ImXp'),1,size(pb,2)).*conj(exp(1j*pb)));
outsum = ImX*vm_kernel(delta,kappa);
outsumW = abs(ImX)*vm_kernel(delta,kappa);
P = (outsum./outsumW).^2;
 
%Von Mises Kernel
function r = vm_kernel(x,kappa)
% C = 1/(2*pi*besseli(0,kappa));
% r = C * exp(kappa*cos(x));
r = exp(kappa*cos(x));