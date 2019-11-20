function [rpd,phasebins] = RPD(x,y,freqVec,width,Fs,numbins)

x = x'; y = y';
Wx = zeros(length(freqVec),size(x,2)); 
Wy = zeros(length(freqVec),size(x,2));

for j=1:length(freqVec)
    Wx(j,:) = energyvec(freqVec(j),detrend(x),Fs,width);
    Wy(j,:) = energyvec(freqVec(j),detrend(y),Fs,width);  
end

X = Wx.*conj(Wy); 
Ph = angle(X);
% Z = zscore(abs(X),[],2);
% iZ = zscore(abs(imag(X)),[],2);
% rZ = zscore(abs(real(X)),[],2);
clear x y Wx Wy X

rpd = zeros(size(Ph,1),numbins);
% rpdhi = zeros(size(Z,1),size(phasebins,2)-1);
% % rpdlo = zeros(size(Z,1),size(phasebins,2)-1);
% rrpdhi = zeros(size(Z,1),size(phasebins,2)-1);
% % rrpdlo = zeros(size(Z,1),size(phasebins,2)-1);
% irpdhi = zeros(size(Z,1),size(phasebins,2)-1);
% irpdlo = zeros(size(Z,1),size(phasebins,2)-1);
% pfr= zeros(size(Z,1),size(phasebins,2)-1);
% rpfr= zeros(size(Z,1),size(phasebins,2)-1);
% ipfr= zeros(size(Z,1),size(phasebins,2)-1);

% for p=1:size(phasebins,2)-1
%     p1 = phasebins(p); p2=phasebins(p+1);
    for f = 1:length(freqVec)
%         pidx = Ph(f,:) >= p1 & Ph(f,:) < p2;
%         rpd(f,p) = nansum(pidx); 
        [rpd(f,:),phasebins] = hist(Ph(f,:),numbins);
%         rpdhi(f,p) = nansum(pidx & Z(f,:)>=2);
% %         rpdlo(f,p) = nansum(pidx & Z(f,:)<=0);
%         rrpdhi(f,p) = nansum(pidx & rZ(f,:)>=2);
% %         rrpdlo(f,p) = nansum(pidx & rZ(f,:)<=0);
%         irpdhi(f,p) = nansum(pidx & iZ(f,:)>=2);
%         irpdlo(f,p) = nansum(pidx & iZ(f,:)<=0);
%         pfr(f,p) = mean(Z(f,pidx));
%         rpfr(f,p) = mean(rZ(f,pidx));
%         ipfr(f,p) = mean(iZ(f,pidx));
    end
% end

% rpd_full = rpd/nansum(nansum(rpd));
rpd = rpd./repmat(nansum(rpd,2),1,size(rpd,2));
% rpdhi = rpdhi./repmat(nansum(rpdhi,2),1,size(rpdhi,2));
% % rpdlo = rpdlo./repmat(nansum(rpdlo,2),1,size(rpdlo,2));
% rrpdhi = rrpdhi./repmat(nansum(rrpdhi,2),1,size(rrpdhi,2));
% % rrpdlo = rrpdlo./repmat(nansum(rrpdlo,2),1,size(rrpdlo,2));
% irpdhi = irpdhi./repmat(nansum(irpdhi,2),1,size(irpdhi,2));
% irpdlo = irpdlo./repmat(nansum(irpdlo,2),1,size(irpdlo,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = energyvec(f,s,Fs,width)
% function y = energyvec(f,s,Fs,width)
%
% Return a vector containing the energy as a
% function of time for frequency f. The energy
% is calculated using Morlet's wavelets. 
% s : signal
% Fs: sampling frequency
% width : width of Morlet wavelet (>= 5 suggested).
%
% 

dt = 1/Fs;
sf = f/width;
st = 1/(2*pi*sf);

t=-3.5*st:dt:3.5*st;
m = morlet(f,t,width);
y = convfft(s,m);
%y = conv(s,m);
%y = (2*abs(y)/Fs).^2;
y = y(ceil(length(m)/2):length(y)-floor(length(m)/2));



function y = morlet(f,t,width)
% function y = morlet(f,t,width)
% 
% Morlet's wavelet for frequency f and time t. 
% The wavelet will be normalized so the total energy is 1.
% width defines the ``width'' of the wavelet. 
% A value >= 5 is suggested.
%
% Ref: Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997)
%
% See also: PHASEGRAM, PHASEVEC, WAVEGRAM, ENERGY 
%
% Ole Jensen, August 1998 

sf = f/width;
st = 1/(2*pi*sf);
A = (st*sqrt(pi))^(-0.5);

y = A*exp(-t.^2/(2*st^2)).*exp(i*2*pi*f.*t);

function g = gaussian(f,t,width)
% function y = morlet(f,t,width)
% 
% Morlet's wavelet for frequency f and time t. 
% The wavelet will be normalized so the total energy is 1.
% width defines the ``width'' of the wavelet. 
% A value >= 5 is suggested.
%
% Ref: Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997)
%
% See also: PHASEGRAM, PHASEVEC, WAVEGRAM, ENERGY 
%
% Ole Jensen, August 1998 

sf = f/width;
st = 1/(2*pi*sf);
A = (st*sqrt(pi))^(-0.5);

g = A*exp(-t.^2/(2*st^2));
