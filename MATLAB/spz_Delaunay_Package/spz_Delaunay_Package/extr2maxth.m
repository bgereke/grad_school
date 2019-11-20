function [x,y] = extr2maxth(M,th)

[C,R] = size(M);

Mid_Mid = zeros(size(M)); % Boolean matrix. True for matrix which min
                          % is at the middle, and max higher than th

for c = 2:C-1   
    for r = 2:R-1
       T = M(c-1:c+1,r-1:r+1) ;
       Mid_Mid(c, r) = (max(max(T)) == T(2, 2)).*(max(max(T))>th);%.*(max(max(T))<th);
    end
end


[x, y] = find(Mid_Mid);
