Tt = 10;
Fs = 2000;
t = (1:Tt*Fs)/Fs;
x = randn(Tt*Fs,1);
nvoices = 32;
freq = 2.^(1:1/nvoices:7);
width = 6;
wx = abs(traces2Wx(x,freq,Fs,width,'Morlet'));
wxsq = wx.^2;
[xZl,yZl] = extr2minth(log(wx), 2);
zt = yZl/Fs;
zf = freq(xZl)';
[T,F] = meshgrid(t,freq);
S = (width+sqrt(2+width^2))./(4*pi*F);
s1 = (width+sqrt(2+width^2))./(4*pi*zf);
E = ones(size(T));
for i=1:length(zt)
%     C=sqrt(2*s1(i)*S./(s1(i)^2+S.^2).*exp(-((T-zt(i)).^2+width^2*(S-s1(i)).^2)./(s1(i)^2+S.^2)));
    C=abs(sqrt(2*s1(i)*S./(s1(i)^2+S.^2)).*exp(-0.5*((T-zt(i)).^2+width^2*(S-s1(i)).^2)./(s1(i)^2+S.^2)+...
        1i*width*(s1(i)+S)./(s1(i)^2+S.^2).*(T-zt(i))));
%     sf = zf(i)/width;st = 1/(2*pi*sf);
%     m = exp(-(t-zt(i)).^2/(2*st^2)).*exp(1i*2*pi*zf(i).*(t-zt(i)));
%     C = abs(traces2Wx(m',freq,Fs,width,'Morlet','area'));
%     C = C/max(C(:));
    E = E.*(1-C).*exp(0);
end

fidx = 105;
[~,idx] = max(wxsq(fidx,2*Fs:Tt*Fs-2*Fs));
tidx = idx+2*Fs;

figure
subplot(4,1,1)
plot(t,wxsq(fidx,:)/wxsq(fidx,tidx),'-k');hold on;plot(t,E(fidx,:)/E(fidx,tidx),'-r');hold off;xlim([2 8])
subplot(4,1,2)
plot(log2(freq),wxsq(:,tidx)/wxsq(fidx,tidx),'-k');hold on;plot(log2(freq),E(:,tidx)/E(fidx,tidx),'-r');hold off;ylim([0 2])
subplot(4,1,3)
surf(t,log2(freq),zeros(size(wx)),wxsq/wxsq(fidx,tidx));view([0 90]);shading interp;axis tight;colormap(gray);colorbar;xlim([2 8]);caxis([0 0.7])
hold on
plot(zt,log2(zf),'.r','MarkerSize',10)
plot(t(tidx)*ones(size(freq)),log2(freq),'-b')
plot(t,log2(freq(fidx))*ones(size(t)),'-b')
subplot(4,1,4)
surf(t,log2(freq),zeros(size(wx)),E/E(fidx,tidx));view([0 90]);shading interp;axis tight;colormap(gray);colorbar;xlim([2 8]);caxis([0 0.7])
hold on
plot(zt,log2(zf),'.r','MarkerSize',10)
plot(t(tidx)*ones(size(freq)),log2(freq),'-b')
plot(t,log2(freq(fidx))*ones(size(t)),'-b')