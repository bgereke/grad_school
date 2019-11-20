
sv0 = [];
sv5 = [];
sv10 = [];

for i=1:12
    [tv,pv,ev,sv] = reconstruct4(0.150,i);

    [pvc,tvc,evc,pfc,tspkc] = collspikes(tv,pv,sv,0);
    s0 =  mean(abs(evc));

    [pvc,tvc,evc,pfc,tspkc] = collspikes(tv,pv,sv,5);
    s5 = mean(abs(evc));

    [pvc,tvc,evc,pfc,tspkc] = collspikes(tv,pv,sv,10); 
    s10 = mean(abs(evc));

    sv0  = [sv0 s0]
    sv5  = [sv5 s5]
    sv10 = [sv10 s10]
end
 
save splitres sv0 sv5 sv10 

