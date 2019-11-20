function [mim_all,mic_all]=cs2nodes_mim(cs,F1,F2)

%function [mim_all,mic_all]=cs2nodes_mim(cs,F1,F2);
%
% Calculates the Multivariate Interaction Measure (MIM)
% and ImCoh with optimzed dipole directions (MIC) for
% all pairs of voxels (nodes inside the brain) on a pre-defined grid
%
% see A. Ewald, L. Marzetti, F. Zappasodi, F. C. Meinecke, and G. Nolte.
% Estimating true brain connectivity from EEG/MEG data invariant to linear
% and static transformations in sensor space. NeuroImage, 60:476 – 488, 2012.
% http://dx.doi.org/10.1016/j.neuroimage.2011.11.084
%
%
% IN   cs        - EEG / MEG sensor level cross-spectrum at a single
%                  frequency bin (might have been averaged over
%                  frequencies); (chans x chans)
%      F1        - inverse filter 1 (chans x voxels x 3)
%      F2        - inverse filter 2 (chans x voxels x 3)
%                  can be used to define reference voxels
%
% OUT  mim_all   - MIM for all voxel pairs (voxels x voxels)
%      mic_all   - maximized ImCoh (MIC) for all voxel pairs
%                  (voxels x voxels)
%
% USAGE
% e.g. to calulate the MIM and the maximized ImCoh (MIC) between all voxel pairs
% defined by the inverse filter F
%      [mim_all,mic_all]=cs2nodes_mim(cs,F,F);
%
% History
% by Guido Nolte - initial version
% updated by Arne Ewald (mail@aewald.net), 12.03.2015 - help inserted, comments added

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

regu=.000001;

if nargin==2
    F2=F1;
end

[nchan, ng, ndipdir]=size(F1);
[~, ng2, ndipdir2]=size(F2);

mim_all=zeros(ng,ng2);
mic_all=zeros(ng,ng2);

csvoxrealinv=zeros(ndipdir,ndipdir,ng);
csvoxrealinvsqrt=zeros(ndipdir,ndipdir,ng);
csvoxvec=zeros(ndipdir,nchan,ng);

for i=1:ng
    Floc=squeeze(F1(:,i,:));
    csloc=Floc'*real(cs)*Floc;
    csvoxrealinv(:,:,i)=inv(csloc+regu*eye(ndipdir)*mean(diag(csloc)));
    csvoxrealinvsqrt(:,:,i)=sqrt(csvoxrealinv(:,:,i));
    csvoxvec(:,:,i)=Floc'*cs;
end


csvoxrealinv2=zeros(ndipdir2,ndipdir2,ng2);
csvoxrealinvsqrt2=zeros(ndipdir2,ndipdir2,ng2);
csvoxvec2=zeros(ndipdir2,nchan,ng2);
for i=1:ng2
    Floc2=squeeze(F2(:,i,:));
    csloc2=Floc2'*real(cs)*Floc2;
    csvoxrealinv2(:,:,i)=inv(csloc2+regu*eye(ndipdir2)*mean(diag(csloc2)));
    csvoxrealinvsqrt2(:,:,i)=sqrtm(csvoxrealinv2(:,:,i));
    csvoxvec2(:,:,i)=Floc2'*cs;
end


for i=1:ng
    %if round(i/100)*100==i;disp(i);end
    csvoxvecloc=csvoxvec(:,:,i);
    caainv=csvoxrealinv(:,:,i);
    caainvsqrt=sqrtm(csvoxrealinv(:,:,i));
    
    for j=1:ng2;
        Floc2=squeeze(F2(:,j,:));
        cab=imag(csvoxvecloc*Floc2);
        cbbinv=csvoxrealinv2(:,:,j);
        X=cab*cbbinv*cab';
        mim_all(i,j)=trace(caainv*X);
        Y=caainvsqrt*X*caainvsqrt;
        [~,s,~]=svd(Y);
        mic_all(i,j)=sqrt(s(1,1));
    end  % j
end % i

