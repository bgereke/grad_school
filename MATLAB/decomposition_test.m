clear all
close all
fclose('all');
Fs = 2000;
win = 1;
t = linspace(0,win,win*Fs);
freq = 2:2:100; width = 7;
X = zeros(length(t),length(t)*length(freq));

%construct dictionary
dt = 1/Fs;
for i = 1:length(freq)
    f = freq(i);
    sf = f/width;
    st = 1/(2*pi*sf);
    tt = linspace(-win/2,win/2,length(t));
    %     A = (st*sqrt(pi))^(-0.5);
    %     y = A*exp(-tt.^2/(2*st^2)).*exp(1i*2*pi*f.*tt);
    y = exp(-tt.^2/(2*st^2)).*exp(1i*2*pi*f.*tt);
    X(:,(i-1)*length(t)+1:i*length(t)) = lagmatrix(y',-length(y)/2:length(y)/2-1);
end
X(isnan(X)) = 0;
X(abs(X)<0.01) = 0;
X = sparse(X);

%contrust test signal
reps = 5;
numwavs = linspace(1,100,10);
mse_m = zeros(length(numwavs),reps);
mse_rr = zeros(length(numwavs),reps);
mse_al = zeros(length(numwavs),reps);
zmse_m = zeros(length(numwavs),reps);
zmse_rr = zeros(length(numwavs),reps);
zmse_al = zeros(length(numwavs),reps);
for r = 1:reps
for n = 1:length(numwavs)
    sidx = randsample(size(X,2),numwavs(n));
    B = zeros(1,size(X,2))';
    B(sidx) = 1;
    % block = floor(length(t)/length(freq));
    % for i = 1:length(freq)
    %     B((i-1)*length(t)+1:i*length(t)) = 1;%/freq(i);
    % %     B((i-1)*length(t)+(i-1)*block+1:(i-1)*length(t)+i*block) = 1;%/freq(i);
    % end
    y = real(sum(X(:,sidx)*B(sidx),2));
%     plot(t,y,'k')
    
    %get ground truth
    gt = zeros(length(freq),length(t));
    for i = 1:length(freq)
        idx = sidx(sidx >= (i-1)*length(t)+1 & sidx <= i*length(t));
        %     idx = (i-1)*length(t)+(i-1)*block+1:(i-1)*length(t)+i*block;
        %     gt(i,:) = abs(sum(X(:,idx)./freq(i),2))';
        gt(i,:) = sum(X(:,idx),2)';
    end
    gt(isnan(gt)) = 0;
    % gt = gt./repmat(max(gt,[],2),1,length(t));
%     figure
    % subplot(2,1,1)
%     imagesc(t,freq,abs(gt));colormap hot;axis xy
    % subplot(2,1,2)
    % plot(freq,mean(gt,2)','-k')
    
    %get morlet convolution
    [mc] = traces2TFR([y y],freq,Fs,width);
    mc(isnan(mc)) = 0;
%     mc = mc./repmat(max(mc,[],2),1,length(t));
%     figure
    % subplot(2,1,1)
%     imagesc(t,freq,abs(mc));colormap hot;axis xy
    % subplot(2,1,2)
    % plot(freq,mean(sqrt(mc),2)','-k')
    % mse_mc(j) = sum(sum((gt-mc).^2));
    
    %get adaptive lasso
    al = zeros(length(freq),length(t));
    rr = zeros(length(freq),length(t));
    spikes = al;
    options.alpha = 0;
    options.nlambda = 15;
    %     tic
    close all
    fclose('all');
    fit1 = glmnet(real(X),y,'gaussian',options);
    options.alpha = 1-eps;
    beta = real(mc./repmat(mean(mc,2),1,length(t)))';
    options.penalty_factor = 1./abs(beta);
    close all
    fclose('all');
    fit2 = glmnet(X,y,'gaussian',options);
    %     toc
    %     figure
    for i = 1:length(freq)
        Brr = zeros(size(X,2),1);
        Brr((i-1)*length(t)+1:i*length(t)) = fit1.beta((i-1)*length(t)+1:i*length(t),end);
%         plot(t,Brr((i-1)*length(t)+1:i*length(t))+i);hold on
        %     al(i,:) = (sum(X*Bal,2))';
        Bal = zeros(size(X,2),1);
        Bal((i-1)*length(t)+1:i*length(t)) = fit2.beta((i-1)*length(t)+1:i*length(t),end);
%         plot(t,Bal((i-1)*length(t)+1:i*length(t))+i);hold on
        al(i,:) = (sum(X*Bal,2))';
        rr(i,:) = (sum(X*Brr,2))';
        spikes(i,:) = Bal((i-1)*length(t)+1:i*length(t));
    end
    mse_m(n,r) = sum(sum((abs(gt)-abs(mc)).^2))/(size(gt,1)*size(gt,2));
    mse_rr(n,r) = sum(sum((abs(gt)-abs(rr)).^2))/(size(gt,1)*size(gt,2));
    mse_al(n,r) = sum(sum((abs(gt)-abs(al)).^2))/(size(gt,1)*size(gt,2));
    zmse_m(n,r) = sum(sum((zscore(abs(gt),[],2)-zscore(abs(mc),[],2)).^2))/(size(gt,1)*size(gt,2));
    zmse_rr(n,r) = sum(sum((zscore(abs(gt),[],2)-zscore(abs(rr),[],2)).^2))/(size(gt,1)*size(gt,2));
    zmse_al(n,r) = sum(sum((zscore(abs(gt),[],2)-zscore(abs(al),[],2)).^2))/(size(gt,1)*size(gt,2));
end
r/reps
end
al(isnan(al)) = 0;
% figure
% al = al./repmat(max(al,[],2),1,length(t));
% subplot(2,1,1)
% imagesc(t,freq,abs(al));colormap hot;axis xy
% subplot(2,1,2)
% plot(freq,mean(al,2)','-k')

figure
plot(numwavs,mean(mse_m,2),'-k');hold on
plot(numwavs,mean(mse_rr,2),'-b');hold on
plot(numwavs,mean(mse_al,2),'-r')

figure
plot(numwavs,mean(zmse_m,2),'-k');hold on
plot(numwavs,mean(zmse_rr,2),'-b');hold on
plot(numwavs,mean(zmse_al,2),'-r')

file1 = 'C:\Data\mouse23\2014-09-29_10-21-24\begin1\CSC5.ncs';
[samples,ts,tt, Fs, bv, ir] = loadEEG2(file1);
Y = bv*samples; clear samples
Y = fftlowpass(Y,Fs,110,120);
tt = tt-min(tt);
[mc] = traces2TFR([Y Y],freq,Fs,width);

B = [];X = [];
dt = 1/Fs;fid = [];
for i = 1:length(freq)
    f = freq(i);
    sf = f/width;
    st = 1/(2*pi*sf);
    y = exp(-(tt-median(tt)).^2/(2*st^2)).*exp(1i*2*pi*f.*(tt-median(tt)));
    posidx = find(mc(i,:) > mean(mc(i,:)));
    lagidx = posidx - round(length(tt)/2);
    x = lagmatrix(y',posidx);
    x = x(:,posidx);
    X = [X sparse(x)];
    fid = [fid freq(i)*ones(1,size(x,2))];
    B = [B mc(i,posidx)];
end
% start = 65;stop = start+1;
% y = Y(start*Fs:stop*Fs-1);
% options.alpha = 0;
% fitr = glmnet(X,y,'gaussian',options);
% options.alpha = 1;
% options.penalty_factor = 1./abs(fitr.beta(:,end));
% %options.cl = [0;Inf];
% tic
% fit = glmnet(X,y,'gaussian',options);
% toc
% % cvglmnetPlot(fit)
% sp = 50;
% al = zeros(length(freq),length(t));
% rr = zeros(length(freq),length(t));
% spikes = al;
% mse = zeros(length(fit.lambda),1);
% for l = 1:length(fit.lambda)
% mse(l) = mean((y-X*fit.beta(:,l)).^2);
% end
% for i = 1:length(freq)
%     Bal = zeros(size(X,2),1);
%     Bal((i-1)*length(t)+1:i*length(t)) = fit.beta((i-1)*length(t)+1:i*length(t),end);
%     plot(tt(start*Fs:stop*Fs-1),Bal((i-1)*length(t)+1:i*length(t))+20*i);hold on
% %     Ber(Ber<1)=0;
%     spikes(i,:) = Bal((i-1)*length(t)+1:i*length(t));
%     al(i,:) = (sum(X*Bal,2))';
%
%     Brr = zeros(size(X,2),1);
%     Brr((i-1)*length(t)+1:i*length(t)) = fitr.beta((i-1)*length(t)+1:i*length(t),end);
%     plot(tt(start*Fs:stop*Fs-1),Brr((i-1)*length(t)+1:i*length(t))+20*i);hold on
%     rr(i,:) = (sum(X*Brr,2))';
% end
% al(isnan(al)) = 0;
% figure
% subplot(2,1,1)
% imagesc(t,freq,abs(al));axis xy;colormap hot%colorbar
% subplot(2,1,2)
% plot(tt(start*Fs:stop*Fs-1),real(sum(al)),'-r')
% hold on;plot(tt(start*Fs:stop*Fs-1),y,'-k')
%
% tic
% [mc] = traces2TFR([y y],freq,Fs,width);
% toc
% figure
% zmc = zscore(mc,[],2);
% zmc(spikes~=0) = -3;
% imagesc(t,freq,zmc,[-2 2]);axis xy
%
% figure
% rrsp = zscore(abs(rr),[],2);
% rrsp(spikes~=0) = -3;
% imagesc(t,freq,rrsp,[-2 2]);axis xy







