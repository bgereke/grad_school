function [C,l] = multihist;  

tetrode = [1 1 1 1 2 2 3 3 3 3 3 3 3 4 4 5 5 6 6 6 6 6 7 7 7 7 7 8 8 8 8 8 8  11 11 11 11 11]; 
cn      = [1 3 6 7 2 3 1 2 4 5 6 7 8 1 5 1 2 3 4 5 6 7 1 2 6 7 8 1 2 3 4 5 10 1  2  4  6  9]; 

lag = 0.2;
dt =  0.001;
bins = 1+2*lag/dt; 
Csum = zeros(1,bins); 
for i=1:length(tetrode)
    for j=1:length(tetrode)
        if i ~= j
            fprintf(1,'%d %d\n',i,j); 
            x1 = loadcell(tetrode(i),cn(i));       
            x2 = loadcell(tetrode(j),cn(j));       
            [C,l] = correlogram2(x1,x2,0.2,0.001);
            Csum = C + Csum; 
        end
    end
end


