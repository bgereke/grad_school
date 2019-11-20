function [Pygx, Px, Pxgy, bins] = tuningcurve(X,y,U)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Pygx, Px, Pxgy] = tuningcurve(X,y,U) returns the nonlinear 
% response of a neuron to the projection of a spatiotemporal stimulus onto
% a single direction in stimulus space.
% X - raw stimulus
% y - neural response
% U - direction in stimulus space (matrix form of STA or eigenvector)
% Pygx - Probability density of a response given a stimulus projection
% Px - Probaility density of a stimulus projection
% Pxgy - Probability density of a stimulus projection given a response
% bins - bins for the stimulus space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

twin = size(U,1);
numpix = size(X,2);

% standardize the stimuli
Xz = (X-repmat(mean(X),size(X,1),1))./repmat(var(X),size(X,1),1);

% project the stimuli onto U

U = U/norm(U);  % normalize direction to unit length
XontoU = zeros(size(X,1)-twin+1,1);

for i=1:size(X,1)-twin+1
    XontoU(i) = reshape(Xz(i:twin+i-1,:),twin*numpix,1)'*...
        reshape(U,twin*numpix,1);
end

% compute the density of XontoU and XontoU given y

numbins = 15;
bins = linspace(min(XontoU),max(XontoU),numbins+1);% create 50 stimulus bins
y = y(twin:end);
Px = zeros(numbins,1);
Pxgy = zeros(numbins,1);

for i = 1:numbins
    Px(i) = length(find(XontoU > bins(i) & XontoU < bins(i+1)));
    Pxgy(i) = sum(y(find(XontoU > bins(i) & XontoU < bins(i+1))));
end

Px = Px/sum(Px);
Pxgy = Pxgy/sum(Pxgy);

% compute the density of y given XontoU

Px(Px==0) = 0.00000000001; %remove 0's from denominator
Pygx = Pxgy./Px;
Pygx = Pygx/sum(Pygx);



