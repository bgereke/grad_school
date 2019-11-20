function w = remove_err(x,maxerr)

w = 0*x;
j = 0;
for i=1:length(x)
    if abs(x(i)) < maxerr 
        j = j + 1;
        w(j) = x(i);
    end
end
w = w(1:j); 
