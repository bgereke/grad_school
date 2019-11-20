function [M] = backfit(y,X,span)

M = zeros(size(X));
alpha = mean(y);
R = y - alpha;
[n,d] = size(X);
fcng = zeros(n,d);
delta = 1;
I = 0;

while delta > 0.001    
    I = I + 1;
%     disp(['Interation number ', int2str(I)]);    
    for i = 1:d        
        J = 1:d;
        rt = R - sum(M(:,J(J~=i)),2);
        fnew = gaussian_kernel(rt,X(:,i),span);
        fnew = fnew-mean(fnew);
        fcng(:,i) = abs(fnew' - M(:,i));
        M(:,i) = fnew;        
    end    
    delta = max(max(fcng));
%     disp(['Max error is ',num2str(delta)]);    
end
    

function [yhat] = gaussian_kernel(y,x,span)

[~,idx] = sort(x);
span = span*length(x);
t = -2.5*span:2.5*span;
k = exp(-t.^2/(2*span^2));
num = convfft(y(idx),k);
num = num(ceil(length(k)/2):length(num)-floor(length(k)/2));
den = convfft(ones(size(y)),k);
den = den(ceil(length(k)/2):length(den)-floor(length(k)/2));
yhat = num./den;
[~,iidx] = sort(idx);
yhat = yhat(iidx);


