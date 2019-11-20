function [c1,c2,c3,c4] = seperate(x)

n = 4;
c1 = [];
c2 = [];
c3 = [];
c4 = [];
 
for i=1:length(x) 
    j = floor(n*x(i,3))+1;
    if j == 1
        c1 = [c1 x(i,2)];
    end
    if j == 2
        c2 = [c2 x(i,2)];
    end
    if j == 3
        c3 = [c3 x(i,2)];
    end
    if j == 4
        c4 = [c4 x(i,2)];
    end
end 
 
