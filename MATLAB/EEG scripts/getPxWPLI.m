function [PxWPLI] = getPxWPLI(x,y,freqVec,width,Fs,abins)

%[PxWPLI,mPxWPLI,zPxWPLI,sPxWPLI] = getPxWPLI(x,y,freqVec,width,Fs,abins)

[xTFR] = traces2TFR([x x],freqVec,Fs,width);
xTFR  =zscore(xTFR,0,2);
% [yTFR] = traces2TFR([y y],freqVec,Fs,width);
% yTFR  =zscore(yTFR,0,2);
[ImX] = traces2ImX(x,y,freqVec,Fs,width);

PxWPLI = zeros(length(abins),length(freqVec));
% mPxWPLI = PxWPLI;
% zPxWPLI = PxWPLI;
% sPxWPLI = PxWPLI;
% nanPFR = PxWPLI;

for a = 1:length(abins)
    PxWPLI(a,:,1) = P_estimator(ImX,xTFR,abins(a),10);
%     PxWPLI(a,:,2) = P_estimator(ImX,yTFR,abins(a),10);
end

% nsh = 100;
% xPxWPLIsh = zeros(size(PxWPLI,1),size(PxWPLI,2),nsh);
% % yPxWPLIsh = zeros(size(PxWPLI,1),size(PxWPLI,2),nsh);
% 
% idx = reshape(randsample(size(ImX,2),size(ImX,2)*nsh,true),nsh,size(ImX,2));
% for s = 1:nsh   
%    for a = 1:length(abins)
%        xPxWPLIsh(a,:,s) =  P_estimator(ImX(:,idx(s,:)),xTFR,abins(a),10);
% %        yPxWPLIsh(a,:,s) =  P_estimator(ImX(:,idx(s,:)),yTFR,abins(a),10);
%    end
%    percent_done = s/nsh*100
% end
% 
% [xci] = prctile(xPxWPLIsh,[2.5 0.5 97.5],3);
% % [yci] = prctile(yPxWPLIsh,[2.5 0.5 97.5],3);
% % nanPFR(:,:,1) = PxWPLI(:,:,1)>xci(:,:,1) & PxWPLI(:,:,1)<xci(:,:,3);
% mPxWPLI(:,:,1) = (PxWPLI(:,:,1) - xci(:,:,2));
% nanPFR(:,:,2) = PxWPLI(:,:,2)>yci(:,:,1) & PxWPLI(:,:,2)<yci(:,:,3);
% mPxWPLI(:,:,2) = (PxWPLI(:,:,2) - yci(:,:,2));
% sPxWPLI = mPxWPLI;
% sPxWPLI(nanPFR) = nan;
% zPxWPLI(:,:,1) = (PxWPLI(:,:,1)-mean(xPxWPLIsh,3))./std(xPxWPLIsh,0,3);
% zPxWPLI(:,:,2) = (PxWPLI(:,:,2)-mean(yPxWPLIsh,3))./std(yPxWPLIsh,0,3);

function P = P_estimator(ImX,TFR,abin,invh)
% edge-corrected kernel density estimator
delta = TFR-abin;
outsum = sum(ImX.*gaussian_kernel(delta,invh),2);
outsumW = sum(abs(ImX).*gaussian_kernel(delta,invh),2);
P = (outsum./outsumW).^2;
 
%Von Mises Kernel
function r = gaussian_kernel(x,invh)
% r =  invh/sqrt(2*pi)*exp(-0.5*x.*x*invh^2);
r =  exp(-0.5*x.*x*invh^2);