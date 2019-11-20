function C2 = morletcorr(TSI,TSJ)

%correlation between wavelets centered at (ti,si) and (tj,sj)
%t is time in seconds
%s is scale
ti = TSI(:,1);
tj = TSJ(:,1);
si = TSI(:,2);
sj = TSJ(:,2);
width = 6;

C2 = abs(sqrt(2*si*sj./(si^2+sj.^2)).*exp(1i*width*(si+sj)./(si^2+sj.^2).*(tj-ti)).*exp(-0.5*((tj-ti).^2+width^2*(sj-si).^2)./(si^2+sj.^2)));