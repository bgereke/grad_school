function [wmcomp,al]=wedgemusic_complete(S,L,vi,P)

%function [wmcomp,al]=wedgemusic_complete(S,L,vi,P)
%
% Performs a Wedge MUSIC scan, i.e. calculates the connectivity (robust to
% volume conduction) between a given reference voxel to all other
% source voxels
%
% A. Ewald, F. S. Avarvand, and G. Nolte., Wedge MUSIC: A novel approach
% to examine experimental differences of brain source connectivity patterns
% from EEG/MEG data. NeuroImage. 101:610–624, ISSN 1095-9572, 2014.
% http://dx.doi.org/10.1016/j.neuroimage.2014.07.011
%
%
% IN   S         - sensor level complex cross-spectrum at a single
%                  frequency bin (might have been averaged over
%                  frequencies); (chans x chans)
%      L         - lead fields (chans x voxels x 3)
%      vi        - voxel index of reference voxels in lead fields (1 x 1)
%      P         - subspace dimension (= estimated number of sources)
%
% OUT  wmcomp    - connectivity of each voxel to reference source
%                  topography v0
%      alpha     - dipole orientations for
%
%
% History
% by Guido Nolte
% updated by Arne Ewald, 12.03.2015 - subspace construction in function
%                                     help inserted, comments added


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

% construct data subspace
[u,s,v]=svd(imag(S));
S=u(:,1:P)*s(1:P,1:P)*v(:,1:P)';


S=S/norm(S,'fro')*sqrt(P);
[~, ngrid, ndim]=size(L);

[u,s,~]=svd(S);
sd=diag(s);
m0=prod(sd(1:P));

Sortho=u(:,1:P);

wmcomp=zeros(ngrid,1);
al=zeros(ndim*2, ngrid);

for i=1:ngrid
    B1=squeeze(L(:,i,:));
    B2=squeeze(L(:,vi,:));
    B1=B1/norm(B1,'fro')*sqrt(ndim);
    B2=B2/norm(B2,'fro')*sqrt(ndim);
    
    nits=10; % to avoid initial value problems in the optimization
    restemp=zeros(nits,1);
    ytemp=zeros(ndim*2, nits);
    m0temp=zeros(nits,1);
    for it=1:nits
        xstart1=randn(3,1);
        xstart2=randn(3,1);
        [restemp(it), ytemp(:,it), m0temp(it)]=estvecvec_wedge(S,B1,B2,P,Sortho,xstart1,xstart2);
    end
    
    [wmcomp(i),idx_its]=min(restemp);
    al(:,i)=ytemp(:,idx_its);
    m0=m0temp(idx_its);
    
end

wmcomp=m0./wmcomp;

