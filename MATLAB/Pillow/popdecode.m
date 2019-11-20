function [x_hat] = popdecode(x,type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x_hat] = popdecode(x,r,'type') implements a Poisson encoding 
% model and decodes a stimulus from the model's population response using 
% the decoding method provided.
% x_hat - stimulus estimate
% x - input stimulus vector (must be between -10 and 10
% type - string specifying the decoding rule. Must be: 'winner-take-all',
% 'center-of-mass', 'template-matching', or 'maximum-likelihood'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Assume encoding model from previous problems
s = -10:0.1:10;
mu = linspace(-10,10,21);
sd = ones(length(mu),1);
a = 15*ones(length(mu),1);
f = zeros(length(mu),length(x));

%Compute expected responses
for i = 1:length(mu)    
    f(i,:) = a(i)*exp(-(x-mu(i)).^2/(2*sd(i)^2));    
end

%Compute Poisson responses
r = poissrnd(f);

%Decode response
if strcmp(type,'winner-take-all')
    [max_r, max_idx] = max(r,[],1);
    x_hat = mu(max_idx);
elseif strcmp(type,'center-of-mass')
    x_hat = mu*r./sum(r);
elseif strcmp(type,'template-matching')
    %Compute templates
    T = zeros(length(mu),length(s));
    for i = 1:length(mu)    
        T(i,:) = a(i)*exp(-(s-mu(i)).^2/(2*sd(i)^2));    
    end
    %Choose x_hat to minimize template sum-squared-errors
    x_hat = zeros(1,length(x));
    for i = 1:length(x)
        E = sum((repmat(r(:,i),1,size(T,2)) - T).^2);
        [min_e, min_idx] = min(E);
        x_hat(i) = s(min_idx);
    end    
elseif strcmp(type,'maximum-likelihood')
    Prgs = zeros(length(s),length(x));
    Ergs = zeros(length(s),length(mu));
    %Compute tuning curves
    for i = 1:length(s)              
        for j = 1:length(mu)    
            Ergs(i,j) = a(j)*exp(-(s(i)-mu(j)).^2/(2*sd(j)^2));
        end
    end
    %Compute respsonse probabilities
    for i = 1:length(s)        
        for j = 1:length(x)
            Prgs(i,j) = prod(exp(-Ergs(i,:)').*(Ergs(i,:)'.^r(:,j))./factorial(r(:,j)));         
        end
    end
    [max_P, max_idx] = max(Prgs,[],1);   
    x_hat = s(max_idx);
end
            
            
   
