function [res,x]=estvec_wedge(A,B,b,P,Aortho,xstart)

% This is a sub-function used for "wedgemusic_scan.m"

% A: data subspace (low rank approximation of ImCs)
% B: Lead field
% b: given reference source topography
% P: number of assumed sources (data subspace size)
% Aortho: eigenvectors spanning data space
% xstart: initital alpha
%
%
% History
% by Guido Nolte
% updated by Arne Ewald

% License
%   Copyright (C) 2015
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

[~, ndipdir]=size(B);

if nargin>4;
    if ~isempty(Aortho);
      u1=Aortho;
    else
      [u1,~,v1]=svd(A);
    end
else
    [u1,~,v1]=svd(A);
end

% no initial value given
if nargin<6    
    xstart=zeros(ndipdir,1);
end

[u2,~,~]=svd([u1(:,1:P),[B,b]]);
O=u2(:,1:P+ndipdir+1);

Asub=O'*A*O;
Bsub=O'*B;
bsub=O'*b;

x=zeros(ndipdir,1);   %#ok<*NASGU>
if nargin>5
  x=xstart;
end

options = optimset('Maxiter', 15, 'Display', 'off', 'Largescale', 'off');
[x, res, ~, outp] = fminunc(@(x) mats2m_fminunc(x, Asub, Bsub, bsub, P), xstart, options);

end

% cost function
function m=mats2m_fminunc(x, A, B, b, P)

Bx=B*x;
C=A+Bx*b'-b*Bx';
[~,s,~]=svd(C);
sd=diag(s);
m=prod(sd(1:P));

end




