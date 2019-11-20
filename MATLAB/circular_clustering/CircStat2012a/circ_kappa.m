function kappa = circ_kappa(alpha,w,dim)
%
% kappa = circ_kappa(alpha,[w])
%   Computes an approximation to the ML estimate of the concentration 
%   parameter kappa of the von Mises distribution.
%
%   Input:
%     alpha   angles in radians OR alpha is length resultant
%     [w      number of incidences in case of binned angle data]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_kappa(alpha, [], dim)
%
%   Output:
%     kappa   estimated value of kappa
%
%   References:
%     Statistical analysis of circular data, Fisher, equation p. 88
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
% Modified by Brian Gereke to handle matrices, 2015

if nargin < 3
  dim = 1;
end

%alpha = alpha(:);

if nargin<2
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
    error('Input dimensions do not match');
  end  
end

N = size(alpha,dim);

if N>1
  R = circ_r(alpha,w,[],dim);
else
  R = alpha;
end

kappa = zeros(size(R));
kappa(R<0.53) = 2*R(R<0.53) + R(R<0.53).^3 + 5*R(R<0.53).^5/6;
kappa(R>=0.53&&R<0.85) = -.4 + 1.39*R(R>=0.53&&R<0.85) + 0.43./(1-R(R>=0.53&&R<0.85));
kappa(R>=0.85) = 1/(R(R>=0.85).^3 - 4*R(R>=0.85).^2 + 3*R(R>=0.85));

if N<15 && N>1
  kappa(kappa<2) = max(kappa(kappa<2)-2*(N*kappa(kappa<2)).^-1,0);    
  kappa(kappa>=2) = (N-1)^3*kappa(kappa>=2)/(N^3+N);
end
