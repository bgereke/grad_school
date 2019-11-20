function [STA,STC,UC,SC,VC,U,S,V] = STA_STC(X,y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [STA,STC,U,S,V] = STA_STC(X,y) returns the "y" triggered
% average stimulus,y-triggered stimulus covariance matrix, along with
% the eigenvectors/values corresponding to the directions in stimulus 
% space causing the greatest excitation/inhibition.
% X - a set of m x n stimuli
% y - trigger variable indexing the stimuli to be included 
% STA - the y-triggered average stimulus
% STC - the y-triggered stimulus covariance matrix
% UC,SC,VC - the singular vectors/values given by svd(C) where C is the 
% covariance matrix taken across all the stimuli
% U,S,V - the singular vectors/values given by svd(STC-C) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numtrigs = sum(y); %number of triggered stimuli
STC_temp = zeros((size(X,1)*size(X,2))^2,1);
STC_n = zeros((size(X,1)*size(X,2))^2,1);
C_temp = zeros((size(X,1)*size(X,2))^2,1);
C_n = zeros((size(X,1)*size(X,2))^2,1);

% standardize the stimuli
Xz = (X-repmat(nanmean(X,3),[1,1,size(X,3)]))./repmat(nanvar(X,[],3),[1,1,size(X,3)]);

% compute STA
trigs(:,:,:) = Xz(:,:,y);
STA = nanmean(trigs,3);

% compute C and STC
for i=1:size(X,3)
    Xc = reshape(Xz(:,:,i),size(X,1)*size(X,2),1);
    C_temp = nansum([C_temp reshape(Xc*Xc',(size(X,1)*size(X,2))^2,1)],2);
    C_n = sum([C_n isnan(reshape(Xc*Xc',(size(X,1)*size(X,2))^2,1))],2);
    if y(i)
        Xstc = reshape(trigs(:,:,i),size(X,1)*size(X,2),1) -...
            reshape(STA,size(X,1)*size(X,2),1);    
        STC_temp = nansum([STC_temp reshape(Xstc*Xstc',(size(X,1)*size(X,2))^2,1)],2);
        STC_n = sum([STC_n isnan(reshape(Xstc*Xstc',(size(X,1)*size(X,2))^2,1))],2);
    end
end

C = reshape(C_temp./C_n,size(X,1)*size(X,2),size(X,1)*size(X,2));
STC = reshape(STC_temp./STC_n,size(X,1)*size(X,2),size(X,1)*size(X,2));

%Do svd
%[UC,SC,VC] = svds(C);
[U,S,V] = svds(STC-C);

