function w = extractphase(S,theta,position)

w = zeros(150,100); 
for i=1:length(S)
    ph = getphase(theta,S(i));
    b  = getbin(position,S(i)); 
    if ph ~= -1
        pl = 1+floor(ph*100); 
        if pl > 100
            pl = 100;
        end
        if pl  < 0
            pl = 1;
        end
        
        w(b,pl) = w(b,pl) + 1;    
    end
end

