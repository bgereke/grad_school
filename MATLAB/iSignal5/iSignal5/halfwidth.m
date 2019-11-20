function [FWHM,ratio]=halfwidth(x,y)
% function [FWHM,ratio]=halfwidth(x,y) computes the approximate full width
% at half maximum (FWHM) of an isolated peak of any smooth shape that has a
% zero baseline. Not highly accurate if the function is too sparsely
% sampled or noisy. Returns the FWHM and the ratio of the lower to upper
% widths, a measure of asymmetry (ratio=1 or -1  for a symmetric peak).
%  Tom O'Haver (toh@umd.edu) April 2016
%
% Example 1:
% x=-5:.1:5;y=sinc(x);
% plot(x,y);
% FWHM=halfwidth(x,y)
%
% Example 2:
% x=[0:.1:10];
% W=3; % W is the true half-width
% y=gaussian(x,5,W);
% FWHM=halfwidth(x,y);
% plot(x,y);
% hold on;plot([5-FWHM/2 5+FWHM/2],[max(y)/2 max(y)/2],'r');hold off
%
format compact
try
    maxy=max(y);
    halfy=maxy/2;
    indmax=round(val2ind(y,maxy));
    n=indmax(1);
    while y(n)>halfy,
        n=n-1;
    end
    y1=interp1([y(n) y(n+1)],[x(n) x(n+1)],max(y)/2);
    n=indmax(1);
    while y(n)>halfy,
        n=n+1;
    end
    y2= interp1([y(n-1) y(n)],[x(n-1) x(n)],max(y)/2);
    FWHM=y2-y1;
    ratio=(y1-x(indmax))./(x(indmax)-y2);
catch
    FWHM=NaN;
end
