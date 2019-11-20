function [x] = mvnRand(mu,COV,n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x] = mvnRand(mu,COV,n) returns n samples from a multivariate
% Gaussian defined by mean mu, covariance COV. 
% mu - mean (must have same number of dimensions as COV
% COV - covariance matrix (must be positive definite)
% n - number of samples
% x - returned samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x1 = normrnd(0,1,n,1);
x2 = normrnd(0,1,n,1);
x = [x1 x2];

[U,S,V] = svd(COV);
A = U*sqrt(S);

x = (A*x'+repmat(mu,1,n))';

