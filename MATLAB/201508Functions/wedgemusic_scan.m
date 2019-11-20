function [wmscan,alpha]=wedgemusic_scan(S,L,v0,P)

%function [wmscan,alpha]=wedgemusic_scan(S,L,v0,P)
%
% Performs a Wedge MUSIC scan, i.e. calculates the connectivity (robust to
% volume conduction) between a given source scalp topography to all other
% source voxels
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
%      L         - lead fields (chans x voxels x 3)
%      v0        - reference source topography (chans x 1)
%      P         - subspace dimension (= estimated number of sources)
%
% OUT  wmscan    - connectivity of each voxel to reference source
%                  topography v0
%      alpha     - dipole orientation at each voxel
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

% normalize data subspace
S=S/norm(S,'fro')*sqrt(P);

% normalize leadfields
[~,ngrid,ndipdir]=size(L); % leadfields: channels x grid points x 3
for i=1:ngrid;
    L(:,i,:)=L(:,i,:)/norm(squeeze(L(:,i,:)),'fro'); %#ok<*AGROW>
end

% normalize reference patterns
v0=v0/norm(v0);

% initial wedge music value (product of singular values)
[u,s,~]=svd(S);
sd=diag(s);
m0=prod(sd(1:P));

Sortho=u(:,1:P);

wmscan=zeros(ngrid,1);
alpha=zeros(ngrid,ndipdir);

for i=1:ngrid
    if ngrid>1;
        B=squeeze(L(:,i,:));
    else
        B=L;
    end
    B=B/norm(B,'fro')*sqrt(ndipdir);
    % Wedge Music at voxel i
    [wmscan(i),y]=estvec_wedge(S,B,v0,P,Sortho);
    alpha(i,:)=y';
end

wmscan=m0./wmscan; % see Eq.(12) in the paper
