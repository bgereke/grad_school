function w  = collspikes_dt(spk,div,bay); 
w = [];
i=0; 
for dt=0.050:0.050:1.500                 
    i = i+ 1;
    phase = 0;
    fname = makename(dt,div,1,bay)  
    load(fname);  
     
    [xvc,pvc,tvc,evcr,pfc,tspkc] = collspikes(xv,tv,pv,sv,spk);
    w = [w mean(abs(evcr))];
end


