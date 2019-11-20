function b = getbin(p,t)

bins = 150;
lgd = 205;    

pos      = getpos(p,t);
x        = pos(1);
y        = pos(2);
l        = getdist([x y]);
b        = round(l*bins/lgd);

if b < 1
    b = 1;
end    
if b > bins
    b = bins;
end


