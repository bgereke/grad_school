function [mim,mic]=cs2reg_mim(cs,F1,F2,p)

%function [mim_all,mic_all]=cs2reg_mim(cs,F1,F2,p);
%
% Calculates the Multivariate Interaction Measure (MIM)
% and ImCoh with optimzed dipole directions for
% all pairs of voxels on a pre-defined grid
%
% IN   cs        - sensor level complex cross-spectrum at a single frequency bin; 
%                  might have been overaged over frequencies
%                 (chans x chans)
%      F1        - inverse filter (chans x voxels x 3) 
%      F2        - inverse filter (chans x voxels x 3) 
%                  can be used to define reference voxels
%      p         - overfitting parameter
%
% OUT  mim_all   - MIM between the two source regions (1 x 1)
%      mic_all   - maximized ImCoh between the two source regions (1 x 1)
%
% USAGE
% e.g. to calulate the MIM and the maximized ImCoh (MIC) 
% between two region with inverse filters F1 & F2
%      [mim,mic]=cs2reg_mim(cs,F1,F2,10);
%
% Arne Ewald, 20.01.2015

regu=.000001;
[nchan, ng, ndipdir]=size(F1);
[~, ng2, ndipdir2]=size(F2);
if nargin<4
warning('Overfitting parameter not set');
p=min(ng,ng2);
end


f1r=reshape(F1, nchan, ndipdir*ng);
f2r=reshape(F2, nchan, ndipdir2*ng2);

% projection of cross-spectra to source space
cs_source=[];
cs_source{1}=f1r'*(cs)*f1r; % aa
cs_source{2}=f1r'*(cs)*f2r; % ab
cs_source{3}=f2r'*(cs)*f2r; % bb
% cs_source{4}=cs{1}';f1r'*(cs)*f1r;

% dimensionality reduction for cross-spectra
cs_red=[];
ca=real(cs_source{1});
[nd,~]=size(ca);
ct=ca+eye(nd,nd)*mean(diag(ca))*10^(-10);  % Regularization
[ua,~,~] = svd(ct);
cs_red{1}=ua(:,1:p)'*cs_source{1}*ua(:,1:p);

cb=real(cs_source{3});
[nd,~]=size(cb);
ct=cb+eye(nd,nd)*mean(diag(cb))*10^(-10);  % Regularization
[ub,~,~] = svd(ct);
cs_red{3}=ub(:,1:p)'*cs_source{3}*ub(:,1:p);

cs_red{2}=ua(:,1:p)'*cs_source{2}*ub(:,1:p);

caainv=inv(real(cs_red{1})+regu*eye(p)*mean(diag(real(cs_red{1}))));
cab=imag(cs_red{2});
cbbinv=inv(real(cs_red{3})+regu*eye(p)*mean(diag(real(cs_red{3}))));
X=cab*cbbinv*cab';
% MIM
mim=(trace(caainv*X));
caainvsqrt=sqrtm(caainv);
Y=caainvsqrt*X*caainvsqrt;
[~,s,~]=svd(Y);
% MIC
mic=sqrt(s(1,1));


end


