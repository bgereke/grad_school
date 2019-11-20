function [STA,STC,U,S,V] = STA_STC(X,y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [STA,STC,U,S,V] = STA_STC(X,y) returns the spike-triggered
% average stimulus,spike-triggered stimulus covariance matrix, along with
% the eigenvectors/values corresponding to the directions in stimulus 
% space causing the great excitation/inhibition.
% X - a T x n stimulus matrix, T (time), n (space)
% y - a length T vector of spike counts
% STA - the spike-triggered average stimulus
% STC - the spike-triggered stimulus covariance matrix
% U,S,V - the singular vectors/values given by svd(STC-C) where C is the 
% covariance matrix taken across all the stimuli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spidx = find(y>0); % spike indices
twin = 10; % num timesteps included in stimulus
spidx = spidx(spidx>twin);
numpix = size(X,2);  % num pixels in stimulus
numspikes = length(spidx);
sptrigs = zeros(twin,numpix,numspikes);  % spike-triggered stimuli
STC = zeros(numpix*twin,numpix*twin);  % spike-triggered covariance matrix
C = zeros(numpix*twin,numpix*twin);  % stimulus covariance matrix


% standardize the stimuli
Xz = (X-repmat(mean(X),size(X,1),1))./repmat(var(X),size(X,1),1);

% compute STA
for i=1:numspikes
    sptrigs(:,:,i) = Xz(spidx(i)-twin+1:spidx(i),:);
end

STA = mean(sptrigs,3);

% compute STC
for i=1:numspikes
    Xi = reshape(Xz(spidx(i)-twin+1:spidx(i),:),numpix*twin,1) -...
        reshape(STA,numpix*twin,1);    
    STC = STC + Xi*Xi';
end

STC = STC/(numspikes-1);

% compute C
for i = twin+1:size(X,1)
    C = C + reshape(Xz(i-twin+1:i,:),twin*numpix,1)*...
        reshape(Xz(i-twin+1:i,:),twin*numpix,1)';
end

C = C/(size(X,1)-twin);
[U,S,V] = svd(STC-C);

