w=0.2*2000;
T = length(x);
W = zeros(T);
x = x/max(x);
% span=0.5/(length(x)/2000);
% x = x-smooth(x,span,'lowess');
for i = 1:T
    idxi = round(i-w/2:i+w/2);
    eidx = idxi>0&idxi<=T;
    idxi = idxi(eidx);
    for j = 1:T
        idxj = idxi+j-i;
        if sum(idxj<=0|idxj>T)==0   
            W(i,j) = exp(-mean((x(idxi)-x(idxj)).^2)/0.02);
        end
    end
end

D = exp(W.^2/10)'*ones(size(x));
smx = exp(W.^2/10)'*x./D;

% x = x/max(x);
% span=0.5/(length(x)/2000);
% x = x-smooth(x,span,'lowess');
% % dx = (x'-x).^2;
% w=0.5*2000;
% win = ones(1,w)/w;
% T = length(x);
% W = zeros(T);
% for i = 1:T-1
%     y = convfft((x(i+1:end)-x(i:end-1)).^2,win);
%     y = y(ceil(length(win)/2):length(y)-floor(length(win)/2));
%     W = W+diag(exp(-y/0.1),i);
% end
% W = W+W';

Ftheta = 8;    % Hz
xbp = fftbandpass(x,Fs,Ftheta-5,Ftheta-4,3*Ftheta,3*Ftheta+1);
tic
T = 4000;
start = 20000;
xidx = start:start+T;
widx = -1000:1000;
weights = exp(-widx.^2/(4*500^2)).*(1-exp(-widx.^2/(2000)));
sampidx = ceil(length(widx)/2):length(x)-ceil(length(widx)/2);
sampidx(sampidx>=start & sampidx<=start+T) = [];
nsamp = 20000;
rsamp = randsample(sampidx,nsamp);
%construct sample patches
spatches = zeros(nsamp,length(widx));
for p = 1:nsamp
   spatches(p,:) = weights.*xbp(widx+rsamp(p)); 
end
%construct x patches
xpatches = zeros(T,length(widx));
for p = 1:T
   xpatches(p,:) = weights.*xbp(widx+xidx(p)); 
end
%compute distances
D = zeros(nsamp,T);
for d = 1:T
    D(:,d) = median((spatches-repmat(xpatches(d,:),size(spatches,1),1)).^2,2);
end
%compute weights
W = exp(-D/0.00025);
Z = ones(1,nsamp)*W;
smx = x(rsamp)'*W./Z;
toc