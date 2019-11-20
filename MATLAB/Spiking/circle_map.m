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

function [map,tmap,cidx] = circle_map(spk_phase,phase,t,kappa,mapAxis)

[map, tmap] = ratemap(spk_phase,phase,t,kappa,mapAxis);
[~, cidx] = (max(map));
% fcen = mapAxis(cidx);

%__________________________________________________________________________
%
% Field functions
%__________________________________________________________________________

% Calculates the rate map.
function [map, tmap] = ratemap(spk_phase,phase,post,kappa,mapAxis)

map = zeros(length(mapAxis),1);
tmap = zeros(length(mapAxis),1);
idx = 0;
for i = mapAxis   
    idx = idx + 1;
    [map(idx), tmap(idx)] = rate_estimator(spk_phase,i,kappa,phase,post);    
end

% Calculate the rate for one position value
function [r, tmap] = rate_estimator(spk_phase,i,kappa,phase,post)
% edge-corrected kernel density estimator
% make (spk_phase-i)->-pi:pi
delta = angle(exp(1j*spk_phase).*conj(exp(1j*i*ones(size(spk_phase)))));
conv_sum = sum(vm_kernel(delta,kappa));
tdiff = diff(post);tdiff(end+1) = median(tdiff);
tdiff(abs(tdiff)>1) = median(tdiff);
delta = angle(exp(1j*phase).*conj(exp(1j*i*ones(size(phase)))));
tmap =  tdiff*vm_kernel(delta,kappa)';
r = conv_sum/tmap;

% Gaussian kernel for the rate calculation
function r = vm_kernel(x,kappa)

%C = 1/(2*pi*besseli(0,kappa));
%r = C * exp(kappa*cos(x));
r = exp(kappa*cos(x));


