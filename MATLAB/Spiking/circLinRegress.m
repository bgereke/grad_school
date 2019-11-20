function [slope, theta_0, circ_corr] = circLinRegress(xx, theta)
%% Based off methods section in Reifenstein et al. PNAS, 2012.
%
%	Inputs - position (xx); phase angle (theta, corrects to radians if necessary)
%
%	Outputs - Regression slope (slope, in radians per unit distance);
%		 Phase intercept (theta_0, in radians);
%		 cicrular linear correlation (circ_corr, range from 0,1);
%

if(max(theta) > 2*pi)
	theta = 2*pi*(theta / 360);
end

ang_bar = @(th) (angle(1/length(th) * sum(exp( 1i * th))));
rr = @(th, ph) sum(sin(th - ang_bar(th)).* sin(ph - ang_bar(ph))) / sqrt(sum(sin(th-ang_bar(th)).^2)*sum(sin(ph-ang_bar(ph)).^2));
R_fit = @(s) -sqrt( (1/length(theta) * sum(cos(theta - s.* xx)))^2 + (1/length(theta) * sum(sin(theta - s.* xx)))^2);

[slope] = fminsearch(R_fit,0);
[theta_0,~,~,~,~] = circularStat(theta(xx <= 0.2*max(xx)));
theta_0=theta_0/(360) * 2 * pi;
% circ_corr = rr(theta, mod(abs(slope)*xx, 2*pi));
circ_corr = rr(theta, xx);