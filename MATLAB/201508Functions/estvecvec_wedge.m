function [res, x, m00]=estvecvec_wedge(A,B1,B2,P,Aortho,xstart1,xstart2)

% This is a sub-function used for "wedgemusic_scan.m"

% A: data subspace (low rank approximation of ImCs)
% B1: Lead field 1
% B2: Lead field 2
% P: number of assumed sources (data subspace size)
% Aortho: eigenvectors spanning data space
% xstart: initital values
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



[~, ndipdir]=size(B1);

if nargin>4;
    if ~isempty(Aortho);
      u1=Aortho;
    else
      [u1,~,~]=svd(A);
    end
else
    [u1,~,~]=svd(A);
end

if nargin<6
    xstart1=zeros(ndipdir,1);
    xstart2=zeros(ndipdir,1);
end

% constructing the low dimensional subspace (p + 2x3)
[u2,~,~]=svd([u1(:,1:P),[B1,B2]]);
O=u2(:,1:P+2*ndipdir);

Asub=O'*A*O;
Bsub1=O'*B1;
Bsub2=O'*B2;

% initializing starting values 
x1=zeros(ndipdir,1);
if nargin>5
  xstart1 =  xstart1/norm(xstart1);
  x1=xstart1;
end
x2=zeros(ndipdir,1);
if nargin>6
  x2=xstart2;
end

m00=mats2m(Asub,Bsub1,Bsub2,0*x1,0*x2,P); % cost function value at 0
m0=mats2m(Asub,Bsub1,Bsub2,x1,x2,P); % cost function value at start values

maxits=1000;

options = optimset('MaxIter', maxits, 'Display', 'off', 'Largescale', 'off');
% Minimization of 5-dim polar coordinates
[xstart1_pol(1), xstart1_pol(2), xstart1_pol(3)] = cart2sph(xstart1(1), xstart1(2), xstart1(3)); 
[x3, res3 , ~, ~, ~, ~] = fminunc(@(x) mats2m_fminunc_pol(x, Asub, Bsub1, Bsub2, P), [xstart1_pol(1:2)'; xstart2], options);

res=res3;

[xx, xy, xz]=sph2cart(x3(1), x3(2), 1);
x=[xx xy xz x3(3:5)']';

end

function m=mats2m(A,B1,B2,x1,x2,P)

Bx1=B1*x1;
Bx2=B2*x2;
C=A+Bx1*Bx2'-Bx2*Bx1';
sd=svd(C);
m=prod(sd(1:P));

end


function m=mats2m_fminunc_pol(x,A,B1,B2,P)

[x1(1), x1(2), x1(3)] = sph2cart(x(1), x(2), 1);
x1 = x1';
x2 = x(3:5);

Bx1=B1*x1;
Bx2=B2*x2;
C=A+Bx1*Bx2'-Bx2*Bx1';
sd=svd(C);
m=prod(sd(1:P));

end



