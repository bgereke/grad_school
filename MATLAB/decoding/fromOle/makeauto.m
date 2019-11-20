function makeauto(x)

lgd = 10; 


ofp = fopen('cgamma.dat','w'); 

for i=1:length(x)-lgd-1
    for j=i:lgd+i
        d = x(j) - x(i); 
        if d < 0.3
            fprintf(ofp,'%d\n',d); 
        end
    end
end
fclose(ofp); 
