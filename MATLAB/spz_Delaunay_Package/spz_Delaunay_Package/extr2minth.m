function [x,y] = extr2minth(M,th)

[R,C] = size(M);

Mid_Mid = zeros(size(M)); % Boolean matrix. True for matrix which min
                          % is at the middle, and max higher than th

for r = 2:R-1   
    for c = 2:C-1
        T = M(max([r-2,1]):min([r+2 R]),max([c-2,1]):min([c+2 C]));
        tr = 3; tc = 3;
        if r == 2
            tr = 2;
        end
        if c == 2
            tc = 2;
        end
        Mid_Mid(r, c) = (min(min(T)) == T(tr,tc)).*(min(min(T))<th);            
    end
end


[x, y] = find(Mid_Mid);

