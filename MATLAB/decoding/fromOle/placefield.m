function w = placefield(tf,p,theta,occ,div)
% function w = placefield(tf,p,theta,occ,div)
%
% Generate the 'place field' for a single cell.
% tf    = [spike-time]  
% p     = [t x y angle]  
% theta = [theta-peaks]
% occ   = [1:bins],  total time at location i

bins = 150;
lgd  = 205; 
w = zeros(bins,div);           
for i=1:length(tf)
    % if mod(i,100) == 0
    %     i/length(tf)
    % end
    phase    = getphase(theta,tf(i)); 
    if phase == 1
        phase = 0.999;
    end
    c        = 1+floor(phase*div);
    pos      = getpos(p,tf(i));
    x        = pos(1);
    y        = pos(2);
    a        = pos(3);   
    l        = getdist([x y]); 
    b        = round(l*bins/lgd); 
          
    if b>0 & b <= bins & phase ~= -1 
        w(b,c) = w(b,c) + 1; 
    end
end 
for i=1:div
    w(:,i) = w(:,i)./occ;
end
