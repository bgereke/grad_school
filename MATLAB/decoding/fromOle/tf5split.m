function [c1,c2,c3,c4,c5] = tf5split(S,theta,position)

c1 = [];
c2 = [];
c3 = [];
c4 = [];
c5 = [];
div = 5; 
for i=1:length(S)
    if mod(i,100) == 0 
        i
    end
    ph = getphase(theta,S(i));
    if ph == 1
        ph = 0.999;
    end
    x  = getbin(position,S(i)); 
    if ph ~= -1
        c  = 1+floor(ph*div);
        if c == 1
            c1 = [c1 [x S(i)]'];
        end
        if c == 2
            c2 = [c2 [x S(i)]'];
        end
        if c == 3
            c3 = [c3 [x S(i)]'];
        end
        if c == 4
            c4 = [c4 [x S(i)]'];
        end
        if c == 5
            c5 = [c5 [x S(i)]'];
        end
    end
end

