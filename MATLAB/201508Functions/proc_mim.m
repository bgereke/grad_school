function [mim,mic,filt_a,filt_b, faall, fball]=proc_mim(cs,na,nb,pa,pb)

% Syntax: 
% function [mim,mic,filt_a,filt_b]=proc_mim(cs,na,nb,pa,pb)
%
% Calculates the Multivariate Interaction Measure (MIM) and the Maximized
% Imaginary Coherency (MIC) between two data spaces A and B. 
%
% MIM can be seen as the total interaction between the spaces A and B
% whereas MIC is the maximized interaction (in terms of the imaginary part
% of coherency).
%
% see A. Ewald, L. Marzetti, F. Zappasodi, F. C. Meinecke, and G. Nolte. 
% Estimating true brain connectivity from EEG/MEG data invariant to linear 
% and static transformations in sensor space. NeuroImage, 60:476 – 488, 2012. 
% http://dx.doi.org/10.1016/j.neuroimage.2011.11.084
%
%
% Input:
%      cs        - complex cross-spectrum for data spaces A and B at a 
%                  single frequency bin (might have been overaged over 
%                  frequencies) having the block structure and size
%
%                                       |-----na-----|------nb-----| 
%                             _         ____________________________
%                             |                      |
%                             na          cs_AA      |    cs_AB
%                             |                      |
%                             _         _____________|______________
%                             |                      |
%                             nb          cs_BA      |    cs_BB
%                             |                      |
%                             _         _____________|______________
%                  
%                   where cs_AA is the cross-spectrum of data space A,
%                   cs_BB the cross-spectrum of data space B, and cs_BB 
%                   the cross-spectrum between data spaces A and B;
%
%                   One will obtain such a structure by combining the data
%                   of the subspaces A and B and calculating the
%                   cross-spectrum for these combined data (see example
%                   below)
% 
%      na        - number of channels in space A
%      nb        - number of channels in space B
%
% Optional input:
%      pa        - overfitting parameter for A
%      pb        - overfitting parameter for B
%
% Output:
%      mim_all   - MIM between the A and B (1 x 1)
%      mic_all   - maximized ImCoh between A & B (1 x 1)
%       filt_a   - spatial filter for space A
%       filt_b   - spatial filter for space B
%
%
% Usage (Example):
%
% nchA=10;  % number of channels in space A
% nchB=15;  % number of channels in space B
% nt=10000; % number of sampling points
% 
% % generate random data
% datA=randn(nchA,nt);
% datB=randn(nchB,nt);
% 
% % combine data in the two spaces
% dat=[datA' datB']';
% 
% % paramters for the estimation of the cross-spectrum (just a toy example)
% segleng=100;segshift=100;epleng=100;maxfreqbin=12;
% 
% % csA=data2cs_event(datA',segleng,segshift,epleng,maxfreqbin);
% % csB=data2cs_event(datB',segleng,segshift,epleng,maxfreqbin);
% 
% % calculation of the cross-spectrum resulting in the appropiate block form
% cs=data2cs_event(dat',segleng,segshift,epleng,maxfreqbin);
% % cs now has size: "(nchA+nchB) X (nchA+nchB) X frequencies"
% 
% % calculation of MIM and MIC at frequency bin #10
% [mim,mic]=proc_mim(cs(:,:,10),nchA, nchB, 5, 7)
%
% History
% Arne Ewald, 11.03.2015: Initial Version

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



regu=.0000000001;
[nall1, nall2]=size(cs);

if (~ismatrix(cs))||(nall1~=nall2)||(isreal(cs))
    error('Wrong format of "cs": must be square, 2D and complex!');
end

if nargin<4
warning('Overfitting parameter not set');
pa=na;
pb=nb;
end

if nall1~=(na+nb)
    error('na+nb must be equal to size of cross-spectrum');
end

cs_aa=cs(1:na,1:na);
cs_ab=cs(1:na,na+1:na+nb);
cs_bb=cs(na+1:end,na+1:end);

% dimensionality reduction for cross-spectrum in space A
ca=real(cs_aa);
[nd,~]=size(ca);
ct=ca+eye(nd,nd)*mean(diag(ca))*10^(-10);  % Regularization
[ua,~,~] = svd(ct);
csr_aa=ua(:,1:pa)'*cs_aa*ua(:,1:pa);


% dimensionality reduction for cross-spectrum in space B
cb=real(cs_bb);
[nd,~]=size(cb);
ct=cb+eye(nd,nd)*mean(diag(cb))*10^(-10);  % Regularization
[ub,~,~] = svd(ct);
csr_bb=ub(:,1:pb)'*cs_bb*ub(:,1:pb);


% dimensionality reduction for cross-spectrum between spaces A & B
csr_ab=ua(:,1:pa)'*cs_ab*ub(:,1:pb);

caainv=inv(real(csr_aa)+regu*eye(pa)*mean(diag(real(csr_aa))));
cab=imag(csr_ab);
cbbinv=inv(real(csr_bb)+regu*eye(pb)*mean(diag(real(csr_bb))));

D=cab*cbbinv*cab'; %#ok<*MINV>

% MIM
mim=(trace(caainv*D));

% [~,s]=eig(caainv*D);
% mim=sum(diag(s));

% spatial filter and MIC
caainvsqrt=sqrtm(caainv);
cbbinvsqrt=sqrtm(cbbinv);
diab=caainvsqrt*cab*cbbinvsqrt;
[u,s,~]=svd(diab*diab');
filt_a=ua(:,1:pa)*caainvsqrt*u(:,1);

faall=ua(:,1:pa)*caainvsqrt*u;

% MIC
mic=sqrt(s(1,1));

[u,~,~]=svd(diab'*diab);
filt_b=ub(:,1:pb)*cbbinvsqrt*u(:,1);

fball=ub(:,1:pb)*cbbinvsqrt*u;



