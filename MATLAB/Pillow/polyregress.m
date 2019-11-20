function [M, w] = polyregress(X,Y,n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [M,w] = polyregress(X,Y,n) finds the coefficients for a degree n 
% polynomial fit to data Y as a function of X
% X = independent variable data points
% Y = dependent variable data points
% n = degree of fit
% M = design matrix 
% w = polynomial coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% construct design matrix

M = zeros(length(X),n+1);

for i = 0:n
    M(:,i+1) = X.^i;
end

% solve regression problem

w = M\Y;
    