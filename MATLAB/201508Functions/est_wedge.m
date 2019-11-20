function [res,x]=est_wedge(A,B,P,Aortho)

% This is a sub-function used for "wedgemusic_scalar.m"

% A: data subspace (low rank approximation of ImCs)
% B: topographies of sources
% P: number of assumed sources (data subspace size)
% Aortho: eigenvectors spanning data space
%
%
% History
% by Guido Nolte
% updated by Arne Ewald, 12.03.2015 - help inserted

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

Bh=B(:,1)*B(:,2)'-B(:,2)*B(:,1)';

if nargin>3
    u1=Aortho;
else
    [u1,~,~]=svd(A);
end

[u2,~,~]=svd([u1(:,1:P),B]);
O=u2(:,1:P+2);

Asub=O'*A*O;
Bhsub=O'*Bh*O;

nite=20;
x=0;
xold=x;
m0=mats2m(Asub,Bhsub,x,P);
m0old=m0;
regu=1;
for i=1:nite;
    [mdiff,mdiffdiff]=mats2malldiff(Asub,Bhsub,x,P);
    dx=-mdiff/(mdiffdiff+regu);
    x=x+dx;
    m0=mats2m(Asub,Bhsub,x,P);
    if m0<m0old;
        regu=regu/10;
        m0old=m0;
        xold=x;
    else
        regu=regu*10;
        x=xold;
    end
end
res=m0old;
return;


function m=mats2m(A,B,x,P)

C=A+x*B;
[~,s,~]=svd(C);
sd=diag(s);
m=prod(sd(1:P));

return;

function [mdiff,mdiff2]=mats2malldiff(A,B,x,P)

dx=.00001;
m3=mats2m(A,B,x,P);
m1=mats2m(A,B,x+dx,P);
m2=mats2m(A,B,x-dx,P);
mdiff=(m1-m2)/(2*dx);
mdiff2=(m1+m2-2*m3)/(dx*dx);

return;


