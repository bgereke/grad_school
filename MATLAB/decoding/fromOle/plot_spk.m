function w = plot_spk;   	  
phase = 1;
cont = 0;
phase = 1;
div = 1;
dt = 0.150;

w = [];
for spk = 0:15                     
    fname = makename(dt,div,phase,cont);  
    load(fname)
    [xvc,pvc,tvc,evc,pfc,tspkc] = collspikes(xv,tv,pv,sv,spk); 
    e1 = mean(abs(evc));  
    w = [w [spk e1 ]']  
end

