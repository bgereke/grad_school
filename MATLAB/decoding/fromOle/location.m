function w = location(p)
% function w = location(p)
%
% Find the probability of the rat being at a given location.
% w : a vector of probabilities
% p : [time x y angle] 
%
% Ole Jensen, May 26, 1998

fprintf(1,'Start: function location \n'); 

bins = 150;
lgd  = 205; 
w = zeros(bins,1);           

for i=1:length(p)
    x        = p(i,2);
    y        = p(i,3);
    a        = p(i,4);   
    l        = getdist([x y]); 
    b        = round(l*bins/lgd); 
          
    if b>0 & b <= bins
        w(b) = w(b) + 1; 
    end
end 

w = w./sum(w);

fprintf(1,'Stop:  function location\n'); 
