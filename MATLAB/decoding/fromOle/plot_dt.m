function w = plot_dt(spk);   	  
phase = 0;
cont = 0;
div = 1;


w = [];
for dt=0.100:0.050:1.500
    fname = makename(dt,div,phase,1)   
    load(fname)
  %  [xvc,pvc,tvc,evc,pfc,tspkc] = collspikes(xv,tv,pv,sv,spk); 
    e1 = mean(abs(ev));  
    w = [w [dt e1 ]']  
end

