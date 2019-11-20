function [maps, pmap, grid, mvl, dbmvl, mvlp, bps, bpt] = ratemap(trans,pos,it,pt,ngrid,kmethod,FWHM,doPerm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute spatial rate maps for calcium flourescence imaging data
%Input:
%trans - [numframes x numcells] matrix of binary calcium transients
%pos - column vector of position samples
%it - column vector of length numframes specifying the time of each imaging frame
%pt - column vector specifying the time of each position sample
%ngrid - number of points at which to evaluate rate map
%kmethod - method for kernel density estimation
%          can be: 'gaussian' or 'vonMises' for 'linear' or 'circular' positions respectively
%FWHM - kernel full width at half max in position units
%doPerm - TRUE/FALSE requests permutation test for p-values and debiasing
%Output:
%maps - [ngrid x numcells] matrix of rate maps
%pmap - ngrid length vector specifying ratemap for transients from entire population
%grid - vector of length ngrid specifying position of rate map
%bps - vector of length numcells specifying spatial information for each map in bits per second
%mvl - mean vector length of transient complex sum (only for vonMises)
%dbmvl - a debiased version of mvl
%mvlp - approximate p-value for each mvl
%bpt - vector of length numcells specifying spatial information for each map in bits per transient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Example - [maps,grid,bps,bpt] = ratemap(trans,lap_position,ct,ard_timestamp,100,'vonMises',30);
%ct and ard_timestamp are from 'plot_dF_on_pos.m'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[numframes, numcells] = size(trans);
mvl = nan(numcells,1);
dbmvl = nan(numcells,1);
mvlp = nan(numcells,1);

%population activity
poptrans = sum(trans,2);

%compute rate maps using specified kernel method
if strcmp(kmethod,'gaussian')
    %set map evaluation grid
    grid = linspace(min(pos),max(pos),ngrid);
    %convert FWHM to gaussian standard deviation
    sd = FWHM/2.355;
    %get position for each imaging frame
    cpos = zeros(size(it));
    for t = 1:length(it)
        [~,midx] = min((pt-it(t)).^2);
        cpos(t) = pos(midx);
    end
    %get rate map for each cell
    delta = repmat(cpos,1,ngrid)-repmat(grid,numframes,1);
    gk = exp(-0.5*delta.*delta/sd^2); %gaussian kernel
    dt = diff(it)';
    dt = [dt(1) dt];
    occupancy = dt*gk;
    maps = (trans'*gk./repmat(occupancy,numcells,1))';
    pmap = (poptrans'*gk./occupancy)';
    
elseif strcmp(kmethod,'vonMises')
    %set map evaluation grid
    grid = linspace(-pi,pi,ngrid);
    %rescale position and FWHM to phase
    pos = pos-(max(pos)-min(pos))/2;
    FWHM = FWHM/max(abs(pos))*pi;
    %convert FWHM to von Mises concentration parameter
    kappa = log(2)/(1-cos(FWHM/2));
    %get position for each imaging frame
    pos = pos/max(abs(pos))*pi;
    cpos = zeros(size(it));
    for t = 1:length(it)
        [~,midx] = min((pt-it(t)).^2);
        cpos(t) = pos(midx);
    end
    %get rate map for each cell
    delta = angle(exp(1i*repmat(cpos,1,ngrid)).*conj(exp(1i*repmat(grid,numframes,1))));
    vmk = exp(kappa*cos(delta)); %von Mises kernel
    dt = diff(it)';
    dt = [dt(1) dt];
    occupancy = dt*vmk;
    px = occupancy/sum(occupancy); %occupancies converted to probabilities
    maps = (trans'*vmk./repmat(occupancy,numcells,1))';
    pmap = (poptrans'*vmk./occupancy)';
    %get mean vector length for each cell
    posrad = exp(1i*cpos);
    [~,midx] = min(abs(repmat(posrad,1,ngrid)-repmat(exp(1i*grid),numframes,1)),[],2);
    weights = 1./px(midx); %inverse occupancies
    for c = 1:numcells
        numt = sum(trans(:,c));
        if numt>3 %don't compute if less than 4 transients
            mvl(c) = abs(weights(trans(:,c)==1)*posrad(trans(:,c)==1)/sum(weights(trans(:,c)==1))); %mean vector length
            if doPerm
                %do permutation testing (num permutations = numframes)
                fmvl  = zeros(numframes,1);
                ftrans = toeplitz([trans(1,c) fliplr(trans(2:end,c)')], trans(:,c)')'; %all temporal shifts of transients
                for f = 1:numframes
                    fmvl(f) = abs(weights(ftrans(:,f)==1)*posrad(ftrans(:,f)==1)/sum(weights(ftrans(:,f)==1)));
                end
                %debias + p-value
                dbmvl(c) = mvl(c) - median(fmvl);
                mvlp(c) = sum(fmvl>=mvl(c))/numframes;
            end
        end
    end
else
    disp('kmethod must be gaussian or vonMises')
end

%get spatial information for each map
dgrid = grid(2) - grid(1);
trate = px*maps; %mean transient rate
bps = dgrid*sum(maps.*log2(maps./repmat(trate,ngrid,1)).*repmat(px',1,numcells)); %bits/sec
bpt = bps./trate; %bits/transient

