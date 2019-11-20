% Calculates the rate map.
function map = ratemap(spkx,spky,posx,posy,post,mapstd,mapgrid)
%spkx - the x locations of the spikes
%spky - the y locations of the spikes
%posx - the x locations of the animal
%posy - the y locations of the animal
%post - the timestamps (in seconds) of the animal position samples
%mapstd - gaussian standard deviation (in position units) for map smoothing
%mapgrid - a vector specifying the x positions at which to evaluate the map
%          assumes square enironment and evaluates at the same positions in
%          the y direction, so positions should be either centered or start
%          at zero. 
%map - length(mapgrid) x length(mapgrid) ratemap evaluated at the locations
%      specified by mapgrid

invh = 1/mapstd;
map = zeros(length(mapgrid),length(mapgrid)); %initilize square map
yy = 0;
%loop over y positions
for mapy = mapgrid
    yy = yy + 1;
    xx = 0;
    %loop over x positions
    for mapx = mapgrid
        xx = xx + 1;
        map(yy,xx) = rate_estimator(spkx,spky,mapx,mapy,invh,posx,posy,post);
    end
end

% Calculate the rate for one position value
function r = rate_estimator(spkx,spky,mapx,mapy,invh,posx,posy,post)

spkcount = sum(gaussian_kernel(spkx-mapx,spky-mapy,invh));
delta = diff(post);delta(end+1) = median(delta);
occupancy =  delta*gaussian_kernel(posx-mapx,posy-mapy,invh)';
r = spkcount/occupancy; 

% Gaussian kernel for the rate calculation
function r = gaussian_kernel(x,y,invh)
% k(u) = ((2*pi)^(-length(u)/2)) * exp(u'*u)
r =  invh^2/sqrt(2*pi)*exp(-0.5*(x.*x*invh^2 + y.*y*invh^2));