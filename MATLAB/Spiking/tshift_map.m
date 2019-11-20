% equalPlot_max('inFile.txt',scale_max)
%
% equalPlot generates ratemaps, plots them and store them to images.
% Multiple session will be read and the maximum rate for a cell will be
% used for plotting in all sessions.
%
% The input file must be on the following format.
%
% C:\Data\TTList.txt
% C:\Data\Begin 1
% C:\Data\Begin 2
% C:\Data\Begin 3
% C:\Data\Begin 4
% and so on ...
%
% The first line specifie the directory to the first session and also the
% name on the t-file list that will be used for all the sessions listed. 
% The t-file list must contain the Mclust t-files, one name for each file. 
% If the list contain cells that only occure in some of the sessions, these
% cells will be plotted as having zero firing rate when missing. Placefield
% images will be stored to both bmp and eps imagefiles to a subdirectory in
% the data folder called placeFieldImages.

function [map,cidx] = tshift_map(aspk_t,invh,tshift)

[map] = ratemap(aspk_t,invh,tshift);
[~, cidx] = (max(map));
% fcen = mapAxis(cidx);

%__________________________________________________________________________
%
% Field functions
%__________________________________________________________________________

% Calculates the rate map.
function [map] = ratemap(aspk_t,invh,tshift)
cut = 3*sqrt(1/invh);
map = zeros(length(tshift),1);
% tmap = zeros(length(tshift),1);
for i = 1:length(tshift)   
    [map(i)] = rate_estimator(aspk_t,tshift(i),invh,cut);    
end

% Calculate the rate for one position value
function [r] = rate_estimator(aspk_t,i,invh,cut)
% edge-corrected kernel density estimator
delta = abs(aspk_t-i);
delta(delta>cut)=[];
% conv_sum = sum(gaussian_kernel(delta,invh));
r = sum(gaussian_kernel(delta,invh));
% tdiff = diff(post);tdiff(end+1) = median(tdiff);
% tdiff(tdiff>1) = median(tdiff);
% delta = abs(atpath-i);
% tdiff(delta>cut) = [];
% delta(delta>cut) = [];
% tmap =  tdiff*gaussian_kernel(delta,invh)';
% r = conv_sum/tmap;

% Gaussian kernel for the rate calculation
function r = gaussian_kernel(x,invh)
% r =  invh/sqrt(2*pi)*exp(-0.5*x.*x*invh^2);
r =  exp(-0.5*x.*x*invh^2);


