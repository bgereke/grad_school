% Example of the RPSEMD performance:
% 
% the algorithm in the manuscript "Regenerated Phase-shifted Sinusoids assisted Empirical Mode Decomposition" by Wang Chenxing, Qian Kemao, Da Feipeng 
% published on IEEE Signal Processing Letters 23(4): 556-560.
%
% -------------------------------------------------------------------------
% Date: Jan 12,2016
% Authors:  W.Chenxing
% For problems with the code, please contact the authors:  
% To:  w.chenxing@gmail.com 
% -------------------------------------------------------------------------
%  This version was run on Matlab R2013b
%--------------------------------------------------------------------------


imf = Rpsemd(Y);  

[cn,co]=size(imf);
% compute OSI
Os=[];
for i=1:cn-1
    Os(i)=sum(imf(i,:).*imf(i+1,:)/co);
end
OSI=sum(Os);

% display
figure; 
subplot(cn+1,1,1); 
plot(Y)
ylabel(['Y']);
set(gca,'xtick',[]);
xlim([1 length(Y)]);
for i=2:cn   
    subplot(cn+1,1,i);  
    plot(imf(i-1,:));  
    ylabel(['c' num2str(i-1)]); 
    set(gca,'xtick',[]);
    xlim([1 length(Y)]);
end;
subplot(cn+1,1,cn+1); 
plot(imf(cn,:));
ylabel(['r']);
xlim([1 length(Y)]);