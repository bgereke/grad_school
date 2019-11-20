function v = makevelo2(p)

c1 = [38 192];
c2 = [132 66];
c3 = [195 209];

v = zeros(length(p)-1,2);  
for i=2:length(p)
    x1 = [p(i-1,2) p(i-1,3)];
    x2 = [p(i,2) p(i,3)];
    a = getarm(x2);
    
    if a == 1
        ac = c2-c1;                                 
    end 
    if a == 2
        ac = c3-c2;
    end 
    if a == 3  
        ac = c1-c3;
    end 
    mc = x2-x1;                      
 
    dproj = ac*mc' * ac / (norm(ac)^2)   ; 
        
    d = 0.4348*norm(dproj); 

    if a == 1 & dproj(1) < 0
        d = -d;
    end
    if a == 2 & dproj(1) < 0
        d = -d;
    end
    if a == 3 & dproj(1) > 0
        d = -d;
    end

    v(i,1:2) = [p(i,1) d/(p(i,1) - p(i-1,1))];
    if mod(i,1000) == 0
        fprintf(1,'%e\n',i/length(p));
    end
end
