function [gim, bias]= proc_gim(cs, p, nave, reg)
%function [gim, bias]= gim(cs, p, nave, reg)
%
% Calculates the Global Interaction Measure
%
% see A. Ewald, L. Marzetti, F. Zappasodi, F. C. Meinecke, and G. Nolte. 
% Estimating true brain connectivity from EEG/MEG data invariant to linear 
% and static transformations in sensor space. NeuroImage, 60:476 – 488, 2012. 
% http://dx.doi.org/10.1016/j.neuroimage.2011.11.084
%
% IN   cs        - cross-spectrum (chans x chans x frequbins)
%      p         - size of subspace prejection (controls overfitting)
%      nave      - (optional) number of averages (e.g. trials or segments)
%                   used to calculate the bias
%      reg       - (optional) regularization parameter (default 0.001)
%
% OUT  gim       - GIM over all frequency bins
%      bias      - bias over all frequency bins (only of "nave" is given)
%
%
% History
% by Arne Ewald (mail@aewald.net), 24.5.2012
% update: regularization added by Arne Ewald, 30.09.2014


% License
%   Copyright (C) 2015  Arne Ewald
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see http://www.gnu.org/licenses/.


if ~exist('nave', 'var') 
    nave=0; 
end

if ~exist('reg', 'var') 
    reg=0.001; 
end

[nch,~, nf]=size(cs);

% fix absolute scaling to avoid nurmerical inconsistencies due to small
% numbers
for f1=1:nf
    cs(:,:,f1)=cs(:,:,f1)./mean(mean(abs(cs(:,:,f1))));
end

uall=cs;

if nargin==1
    warning('Overfitting parameter not set');
    p=nch;
end

% singular value decomposition for dimensionality reduction (overfittig)
for f1=1:nf;
    ci=imag(cs(:,:,f1));
    cit=ci*ci'+eye(nch,nch)*mean(diag(ci*ci'))*10^(-10);  % Regularization
    [u,~,~] = svd(cit);
    uall(:,:,f1)=u;
end

gim=zeros(nf,1);

for f1=1:nf;
     csloc=uall(:,1:p,f1)'*cs(:,:,f1)*uall(:,1:p,f1); % Dimensionality reduction
     csloc=csloc+eye(p,p)*(mean(diag(real(csloc)))*reg); % Regularization
     csloc = cs(:,:,f1);
%      keyboard
     gim(f1)=0.5*trace(inv(real(csloc))*imag(csloc)*...
                                   inv(real(csloc))*(imag(csloc))');
%      e = eig(inv(real(csloc))*imag(csloc)*...
%                                    inv(real(csloc))*(imag(csloc))');
%      gim(f1) = max(e);
end % for f


% Calculate Bias
bias=zeros(1,1);
if nave~=0
    bias=(p*(p-1))/(4*nave);    
end
