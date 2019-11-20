function [d] = mvnDensity(x,mu,COV)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [d] = mvnDensity(x,mu,cov) returns the value of the density of a
% multivariate Gaussian, defined by mean x and covariance matrix COV,
% evaluated at the point x.
% x - evaluation point (must be the same dimension as COV)
% mu - means of the Guassian along the dimensions of x
% COV - covariance matrix of the Gaussian
% d - value of the density evaluated at x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = length(x);
d = 1/sqrt((2*pi)^k*det(COV))*...
    exp(-0.5*(x-mu)'*COV^-1*(x-mu));

