function  nall=recon(fall,p)
%
%  function  nall=recon(fall,p)
%




%------------------------------------------------------------
% Define constants
%------------------------------------------------------------
bins = 150; 
tmin = 1e4;               
tmax = 1.1e4;


%------------------------------------------------------------
% Truncate data files (tmin to tmax)           
%------------------------------------------------------------


maxind = 0;
minind = 1e4;
for i=1:size(fall,2) 
     ind = find(fall(:,i) >= tmin &  fall(:,i) <= tmax);
     if min(ind) < minind
         minind = min(ind);
     end  
     if max(ind) > maxind
         maxind = max(ind);
     end  
end

fall = fall(minind:maxind,:);

ind = find(p(:,1) >=tmin &  p(:,1) <= tmin + (tmax-tmin)/2);
p = p(ind,:);


nall = zeros(bins,size(fall,2));
for j=1:size(fall,2) 
    fprintf('%d ',j);
    for i=1:size(fall,1) 
        ts = fall(i,j); 
        if ts > tmin & ts < (tmin + (tmax-tmin)/2)
            xt  = getbin(p,fall(i,j)); 
            nall(xt,j) = nall(xt,j) + 1; 
        end
    end
end

nall = filtplacefields(nall); 
