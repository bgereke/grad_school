function [A,M] = Add_For_Reg(y,X,k,span,e)

if size(y,1) == 1
    y = y';
end
if size(X,1) ~= size(y,1)
    X = X';
end

A = zeros(1,size(X,2));
D  = 1:size(X,2);
M = zeros(size(X));
alpha = mean(y);
R = y - alpha;

for kk = 1:k
    disp(['Interation number ', int2str(kk)]);
    for d = D
        M(:,d) = gaussian_kernel(R,X(:,d),span);
    end    
    [~,midx] = max(abs(M(:,D)'*R));
    A(kk) = D(midx);
    D(midx) = [];
    M(:,sort(A(1:kk))) = backfit(y,X(:,sort(A(1:kk))),span);
    rt = y - alpha - sum(M,2);
    dMISE = norm(R)-norm(rt);
    R = rt;
    disp(['dMISE is ',num2str(dMISE)]);
%     if dMISE <= e
% %         break
%     else
%         R = rt;
%     end
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