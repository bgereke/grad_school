function z = thetacorr(c,j)
z =zeros(1,100);

z(1) = 0; 
k = 1;
for i=1:100  %length(c)
    if c(i,3) == j
        z(k) = c(i,2);
        k= k + 1;
    end
end
