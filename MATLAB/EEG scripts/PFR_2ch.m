function [medPFRznorm] = PFR(S,Sp,freqVec,width,Fs,phasebins)

[theta_phase, thetaBP] = thetaphase(Sp,Fs);
[TFR,timeVec,freqVec] = traces2TFR([S S],freqVec,Fs,width);
% [TFRInfnorm, TFRcounts] = histTFR(TFR,freqVec,0:0.1:1);
TFRznorm = zscore(TFR,0,2);
% TFR2norm = norm2(TFR);
%TFRperc = zeros(size(TFR));

% for i=1:size(TFR,1)
%     TFRperc(i,:) = tiedrank(TFR(i,:))/size(TFR,2);
% end

% medPFR2norm = zeros(size(TFR,1),size(phasebins,2)-1);
%medPFRznorm = zeros(size(TFR,1),size(phasebins,2)-1);
% medPFRInfnorm = zeros(size(TFR,1),size(phasebins,2)-1);
medPFRperc = zeros(size(TFR,1),size(phasebins,2)-1);
%varPFR = medPFR;
percbins = 0:0.05:1;
for p=1:size(phasebins,2)-1
    
    pidx = theta_phase >= phasebins(p) & theta_phase < phasebins(p+1);
    
%     medPFRperc(:,p) = sum(TFRperc(:,pidx) > 0.95,2)/size(pidx,2);
    
%     medPFR2norm(:,p) = mean(TFR2norm(:,pidx),2);
    medPFRznorm(:,p) = mean(TFRznorm(:,pidx),2);
%     medPFRInfnorm(:,p) = mean(TFRInfnorm(:,pidx),2);
%     medPFRperc(:,p) = median(TFRperc(:,pidx),2);
%     varPFR(:,p) = var(TFRnorm(:,pidx),0,2);
%     stdPFR(:,p) = std(TFRnorm(:,pidx),0,2);

%     if p == 10
%         figure
%         hist(TFRperc(5,pidx),percbins);
%     elseif p == 35
%         figure
%         hist(TFRperc(5,pidx),percbins);
%     end
    
end

percbins = 0:0.1:1;

% figure
% imagesc(phasebins(1:end-1),freqVec(1:end),medPFR2norm(1:end,:)); axis xy; colorbar
% xlabel('phase (rads)');ylabel('frequency (Hz)');title('2-norm');
% figure
% imagesc(phasebins(1:end-1),freqVec(1:end),medPFRznorm(1:end,:)); axis xy; colorbar
% xlabel('phase (rads)');ylabel('frequency (Hz)');title('z-norm');
% figure
% imagesc(phasebins(1:end-1),freqVec(1:end),medPFRInfnorm(1:end,:)); axis xy; colorbar
% xlabel('phase (rads)');ylabel('frequency (Hz)');title('Inf-norm');
% figure
% imagesc(phasebins(1:end-1),freqVec(1:end),medPFRperc(1:end,:)); axis xy; colorbar
% xlabel('phase (rads)');ylabel('frequency (Hz)');title('Percentiles');
% figure
% imagesc(timeVec,freqVec,TFR2norm);axis xy;xlabel('time');ylabel('frequency (Hz)');title('2-norm');
% figure
% imagesc(timeVec,freqVec,TFRznorm);axis xy;xlabel('time');ylabel('frequency (Hz)');title('z-norm');
% figure
% imagesc(timeVec,freqVec,TFRInfnorm);axis xy;xlabel('time');ylabel('frequency (Hz)');title('Inf-norm');


function y = norm2(x)

y = zeros(size(x));
norms = zeros(size(x,1));

for i=1:size(x,1)
    
    norms(i) = norm(x(i,:));
    y(i,:) = x(i,:)/norms(i);
    
end

y;

    