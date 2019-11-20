% DETECT_EXTREMA: detect the indexes of extrema, including minima and maxima
%
% INPUT:
% Y: one-dimensional data
%
% OUTPUT:
% indmax: indexes of maxima
% indmin: indexes of minima
% indminmax: indexes of all extrema
% 
% Refer to emd.m written by G. Rilling: http://perso.ens-lyon.fr/patrick.flandrin/emd.html 

function [indmax,indmin,indminmax]=Detect_extrema(Y)

nY = length(Y);

% the difference of the data
dY = diff(Y);
if dY==0
    indmax=[];indmin=[];indminmax=[];return;
end

ndY = length(dY);
dY1 = dY(1:ndY-1);     dY2 = dY(2:ndY);

% the indexes of minima and maxima
indmin = find(dY1.*dY2<0 & dY1<0)+1;   
indmax = find(dY1.*dY2<0 & dY1>0)+1;

%% refine if adiacent values are the same
if any(dY==0)
  
  imax = [];  imin = [];  r_flat1=[];  l_flat1=[];
  
  % find the local area like a basin
  flat = (dY==0);             % =1, the same values
  dflat = diff([0 flat 0]);     % enlarge the length to match the length of Y    
  l_flat = find(dflat == 1);       % the left border of the basin
  r_flat = find(dflat == -1);      % the right border of the basin
  
  % exclude the left border if it is also the left border of a basin area
  if l_flat(1) == 1   % left border of a basin area = left border of Y ?
      if length(l_flat) == 1 
          l_flat = []; r_flat1 = r_flat(1); r_flat = [];
      else
          l_flat = l_flat(2:end); r_flat1 = r_flat(1); r_flat = r_flat(2:end);
      end
  end
   
  % exclude the right border if it is also the right border of a basin area
  if length(l_flat) > 0    % exist a whole basin area
      if r_flat(end) == nY  % right border of basin area = right border of Y ?
        if length(l_flat) == 1
            l_flat1 = l_flat(1); l_flat = [];  r_flat = [];
        else
            l_flat1 = l_flat(1); l_flat = l_flat(1:(end-1));  r_flat = r_flat(1:(end-1));                    
        end
      end
  end
  
  nl_flat = length(l_flat);
  % refinement
  if nl_flat > 0
    for k = 1:nl_flat
      if dY(l_flat(k)-1) > 0   % the flat peak ?
        if dY(r_flat(k)) < 0
          imax = [imax round((r_flat(k)+l_flat(k))/2)];
        end
      else                % the flat valley ?
        if dY(r_flat(k)) > 0
          imin = [imin round((r_flat(k)+l_flat(k))/2)];
        end
      end
    end
  end
  if ~isempty(r_flat1)
      if dY(r_flat1)<0
          imax = [r_flat1, imax];
      else
          imin = [r_flat1, imin];
      end
  end
  if ~isempty(l_flat1)
      if dY(l_flat1-1)<0
          imin = [imin, l_flat1];
      else
          imax = [imax, l_flat1];
      end
  end  
  
  % range the new maxima and minima
  if length(imax) > 0
    indmax = sort([indmax imax]);
  end
  if length(imin) > 0
    indmin = sort([indmin imin]);
  end  
end  
%%
% mix the extrema
indminmax=sort([indmin indmax]);
end