function [Cxy,lag] = correlogram(x1,x2,maxlag,binsize);
% function [Cxy,lag] = correlogram(x1,x2,maxlag,binsize);
%
% Calculate the crosscorelogram, where x1 and x2 denote vectors
% of subsequent discrete events. The routine is optimized for
% autocorrelograms. 

lag = (-floor(maxlag/binsize) : floor(maxlag/binsize))*binsize; 
Cxy = zeros(size(lag)); 
% h = waitbar(0,'Please wait...');

j0 = 1;
for i=1:length(x1)

    if (mod(i,100) == 0) 
%        waitbar(i/length(x1))   
    end
    while j0 < length(x2) & x2(j0) < x1(i)
        j0 = j0 + 1; 
    end
    j = j0;  
    while j < length(x2) & abs(x2(j) - x1(i)) < maxlag
        l = x2(j) - x1(i);
        p = floor((maxlag+l)/binsize)+1;  
        Cxy(p) = Cxy(p) + 1;
        j = j + 1;
    end
    j = j0-1;
    while j > 0 & abs(x2(j) - x1(i)) < maxlag
        l = x2(j) - x1(i);
        p = floor((maxlag+l)/binsize)+1;  
        Cxy(p) = Cxy(p) + 1;
        j = j - 1;
    end
end
% close(h)
Cxy(length(Cxy)) = Cxy(length(Cxy)-1); 
if (length(x1) == length(x2)) Cxy(floor(length(Cxy)/2)+1) = 0; end;

if (nargout == 0)
   plot(lag,Cxy) , grid on
   xlabel('lag'), ylabel('correlation');
end



