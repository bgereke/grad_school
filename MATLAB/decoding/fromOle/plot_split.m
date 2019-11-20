function w = plot_split(spk,bay);   	  
phase = 1;
cont = 0;
phase = 1;
div = 1;
dt = 0.150;

w = [];
for div=1:9                         
    fname = makename(dt,div,1,bay);  
    load(fname)
    [xvc,pvc,tvc,evc,pfc,tspkc] = collspikes(xv,tv,pv,sv,spk); 
    e1 = mean(abs(evc));  
    w = [w [div e1 ]']  
end

