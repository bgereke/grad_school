function w = tmp(S,theta)

for i=1:length(S)
    if mod(i,100) == 0 
        i
    end
    ph = getphase(theta,S(i));
    if ph ~= -1
        w = [w ph]; 
    end
end

