function [hout wout] = virmenArrangeGrid(h,w,n)

wout = 1;
hout = 1;
while n > wout*hout
    if h/(hout+1) < w/(wout+1)
        wout = wout+1;
    else
        hout = hout+1;
    end
end