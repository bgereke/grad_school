function [quadt,quadf,quadw,quady] = heisenbergquad(time,freq,roots,stdev,width)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%generates quadrature points corresaponding to a heisenberg tiling of the
%time-frequency plane for a Morlet cwt 
%
%Inputs:
%time - time vector in seconds at which the cwt is evaluated
%freq - frequency vector in Hz at which the cwt is evaluated
%roots - 2 x n matrix of cwt roots, top row is time in sec, 
%        bottom row is frequency in Hz 
%stdev - quadrature spacing in time/frequency resolution deviations
%        (assumed equal for time and frequency)
%width - width parameter used in Morlet cwt
%
%Outputs:
%quadt - vector of quadrature time points
%quadf - vector of quadrature frequency points
%quadw - quadrature weight associated with each point
%        (area of tile)/(# of points in tile)
%quady - quadrature pseudodata I(zero)./quadw where I(zero) is the
%        indicator function specifying if the point is a cwt zero or not
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%rescale frequency into std's
fscaled = cumsum(gradient(freq)*2*width./freq);

%get quadrature frequencies
fquads = pchip(fscaled,freq,ceil(min(fscaled)):stdev:floor(max(fscaled)));

%get maximum time deviations for each quadrature frequency
maxt = max(time);
quadt_max = maxt/width*2*pi.*fquads;

%get the quadrature points
quadt = []; quadf = []; quada = [];
for f=1:length(fquads)
    %quadrature times
    tquads = 0.5*stdev:stdev:quadt_max(f);
    
    %quadrature points
    quadt = [quadt tquads*width./(2*pi*fquads(f))];
    quadf = [quadf fquads(f)*ones(size(tquads))];
    
    %quadrature areas (sd_t*sd_f = 1/(4*pi) for Morlet cwt)
    if f ==1
        %area for lower border tiles
        quada = [quada stdev*(stdev/2+ceil(min(fscaled))-min(fscaled))/(4*pi)*ones(size(tquads))];
        %area for lower right border tile
        quada(end) = (stdev/2+ceil(min(fscaled))-min(fscaled))*(stdev/2+quadt_max(f)-max(tquads))/(4*pi);
    elseif f == length(quadf)
        %area for upper border tiles
        quada = [quada stdev*(stdev/2+max(fscaled)-floor(max(fscaled)))/(4*pi)*ones(size(tquads))];
        %area for upper right border tile
        quada(end) = (stdev/2+max(fscaled)-floor(max(fscaled)))*(stdev/2+quadt_max(f)-max(tquads))/(4*pi);
    else
        %area for the full tiles 
        quada = [quada stdev^2/(4*pi)*ones(size(tquads))];
        %area for the right border tiles
        quada(end) = stdev*(stdev/2+quadt_max(f)-max(tquads))/(4*pi);
    end
end

%determine if any quadrature points and roots are exactly equal
roott = roots(1,:);rootf = roots(2,:);
Iquad = zeros(size(quadt)); %indicator function for the quadrature points
Iroots = zeros(size(roott));
for r = 1:length(roott)
    Iqt = roott(r)==quadt;
    Iqf = rootf(r)==quadf;
    if any(Iqt) && any(Iqf)
        Iquad(Iqt && Iqf) = 1;
        Iroots(r) = 1;
    end
end
roott(logical(Iroots)) = [];
rootf(logical(Iroots)) = [];

%get the # of quadrature points + roots in each tile
nearest_f = zeros(size(roott));
nearest_t = zeros(size(roott));
for r = 1:length(roott)
   [~,fidx] = min(abs(log2(rootf(r))-log2(quadf)));
   nearest_f(r) = quadf(fidx);
   quadt_f = quadt(quadf==nearest_f(r));
   [~,tidx] = min(abs(roott(r)-quadt_f));
   nearest_t(r) = quadt_f(tidx);
end

numquad = ones(size(quadt));
for q = 1:length(quadt)
   numquad(q) = 1 + sum(nearest_t==quadt(q) & nearest_f==quadf(q)); 
end
numroot = ones(size(roott));roota = ones(size(roott));
for r = 1:length(roott)
   numroot(r) = 1 + sum(nearest_t==nearest_t(r) & nearest_f==nearest_f(r)); 
   roota(r) = quada(quadt==nearest_t(r) & quadf==nearest_f(r));
end

quadt = [roott quadt];
quadf = [rootf quadf];
quadw = [roota./numroot quada./numquad];
quady = [ones(size(roott)) Iquad]./quadw;

