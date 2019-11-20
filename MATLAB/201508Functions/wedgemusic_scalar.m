function [wmscalar,alpha]=wedgemusic_scalar(S,t1,t2,P)

%function [wmscalar,alpha]=wedgemusic_scalar(S,t1,t2,P)
%
% Performs scalar Wedge MUSIC, i.e. calculates the connectivity (robust to
% volume conduction) between two given source scalp topographies t1 and t2
%
% A. Ewald, F. S. Avarvand, and G. Nolte., Wedge MUSIC: A novel approach
% to examine experimental differences of brain source connectivity patterns
% from EEG/MEG data. NeuroImage. 101:610–624, ISSN 1095-9572, 2014.
% http://dx.doi.org/10.1016/j.neuroimage.2014.07.011
%
%
% IN   S         - sensor level complex cross-spectrum at a single%
%                  frequency bin (might have been averaged over
%                  frequencies); (chans x chans)
%      t1        - first source topography(chans x 1)
%      t2        - second source topography(chans x 1)
%      P         - subspace dimension (= estimated number of sources)
%
% OUT  wmscalar    - connectivity of each voxel to reference source
%                  topography t1
%      alpha     - dipole orientation at each voxel
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

xfac=norm(S,'fro')/sqrt(P);

S=S/norm(S,'fro')*sqrt(P);
[~, ngrid]=size(t2);

t1=t1/norm(t1);

[~,s,~]=svd(S);
sd=diag(s);
m0=prod(sd(1:P));


wmscalar=zeros(ngrid,1);
alpha=zeros(ngrid,1);
for i=1:ngrid
    loctopo=t2(:,i);
    loctopo=loctopo/norm(loctopo);
    B=[t1,loctopo];
    [wmscalar(i),y]=est_wedge(S,B,P);
    alpha(i)=y;
end

for i=1:ngrid;
    alpha(i)=alpha(i)*xfac/(norm(t2(:,i))*norm(t1));
end

wmscalar=m0./wmscalar;
