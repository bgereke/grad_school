function [MChist] = morletcorrhist(rt,rf,t,f,fgrid,width)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Generates the Morlet correlation function C(t,f) induced by cwt roots
%observed at all times t_0:t and frequencies leading up to time t 
%
%Inputs:
%rt - time vector in seconds at which the cwt roots are observed
%rf - frequency vector in Hz at which the cwt roots are observed
%t - time vector in seconds at which data points are taken 
%f - frequency vector in Hz at which data points are taken
%fgrid - frequency grid in Hz at which to evaluate correlation function
%        for each of the time points in t
%width - width parameter used in Morlet cwt
%
%Outputs:
%MChist - length(fgrid) x length(t) matrix giving the history induced 
%         correlation function from points (rt, rf) for each point (t, f)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set vector orientations
if size(t,1) > size(t,2)
    t = t';
end
if size(fgrid,1) < size(fgrid,2)
    fgrid = fgrid';
end

%convert frequency to scale
s = (width+sqrt(2+width^2))./(4*pi*fgrid);
rs = (width+sqrt(2+width^2))./(4*pi*rf);

%compute the correlation function
MChist = ones(length(fgrid),length(t));
for r=1:length(rt) %loop through each of the roots
    %set logical indices s.t. roots influence current and future time points, 
    %but not themselves
    tidx = t>=rt(r) & f~=rf(r);
    MChist(:,tidx) = MChist(:,tidx).*...
        (1-abs(sqrt(2*rs(r)*s./(rs(r)^2+s.^2)).*...
        exp(1i*width*(rs(r)+s)./(rs(r)^2+s.^2).*...
        (t(tidx)-rt(r))).*exp(-0.5*((t(tidx)-rt(r)).^2+width^2*...
        (s-rs(r)).^2)./(rs(r)^2+s.^2))));
end