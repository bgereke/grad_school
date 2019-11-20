function err = geterror(xp,xt)
lgd = 205; 
bins = 150; 

    if xp > xt
        l1 = xp - xt;
        l2 = l1 - bins;
    else
        l1 = bins - (xt - xp);
        l2 = l1 - bins;
    end
    if abs(l1) < abs(l2)
        err = l1;
    else
        err = l2;
    end


err = 205*err/150; 
