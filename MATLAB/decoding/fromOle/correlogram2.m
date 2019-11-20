function [Cxy,lag] = correlogram(x1,x2,maxlag,binsize);
% function [Cxy,lag] = correlogram(x1,x2,maxlag,binsize);
%
% Calculate the crosscorelogram, where x1 and x2 denote vectors
% of subsequent discrete events. The routine is optimized for
% autocorrelograms. 

lag = (-floor(maxlag/binsize) : floor(maxlag/binsize))*binsize; 
Cxy = zeros(size(lag)); 
% h = waitbar(0,'Please wait...');
bins = length(lag);                          


j0 = 1;
for i=1:length(x1)
    while j0 < length(x2) &  x2(j0) < x1(i)  
        j0 = j0 + 1; 
    end
    j = j0;  
    while j <= length(x2) & abs(x2(j) - x1(i)) < maxlag
        l = x2(j) - x1(i);
        p = 1+floor(bins*(maxlag + l)/(2*maxlag)); 
        if p <= bins & p > 0
            Cxy(p) = Cxy(p) + 1;
        end
        j = j + 1;
    end
    j = j0-1;
    while j > 0 & abs(x2(j) - x1(i)) < maxlag
        l = x2(j) - x1(i);
        p = 1+floor(bins*(maxlag + l)/(2*maxlag)); 
        if p <= bins & p > 0
            Cxy(p) = Cxy(p) + 1;
        end
        j = j - 1;
    end
end
% close(h)
if length(x1) == length(x2) 
    if x1 == x2
        Cxy(floor(length(Cxy)/2)+1) = 0; 
    end
end;

if (nargout == 0)
   plot(lag,Cxy) , grid on
   xlabel('lag'), ylabel('correlation');
end



