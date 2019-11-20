numrewards = 6; %number of rewards 
prop = 0.5; %von mises width expressed as proportion of track
numpos = 384; %number of discrete positions
center = round(0.75*numpos); %center of reward location
FWHM = prop*numpos; %full width at half max
pos = 1:numpos; %positions in squares
radpos = linspace(-pi,pi,numpos); %positions in radians
%rescale position and FWHM to phase
pos = (pos-(max(pos)-min(pos))/2);
FWHM = FWHM/max(abs(pos))*pi;
%convert FWHM to von Mises concentration parameter
kappa = log(2)/(1-cos(FWHM/2));
%make von mises distribution
delta = angle(exp(1i*repmat(radpos(center),1,numpos)).*conj(exp(1i*radpos)));
vmk = exp(kappa*cos(delta));
vmk = vmk/sum(vmk); %make probability distribution
%generate a reward sample for a single lap
rewsamp = randsample(1:numpos,numrewards,true,vmk);