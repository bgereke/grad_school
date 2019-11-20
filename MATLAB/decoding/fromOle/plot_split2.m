function [er,ep] = plot_split(spk);   	  
dt = 0.150;
bay = 1;

er = [];
ep = [];
for div=2:9                         
    phase = 0; 
    fname = makename(dt,div,phase,bay);  
    load(fname)
    [xvc,pvc,tvc,evc,pfc,tspkc] = collspikes(xv,tv,pv,sv,spk); 
    er = [er mean(abs(evc))];  

    phase = 1; 
    fname = makename(dt,div,phase,bay);  
    load(fname)
    [xvc,pvc,tvc,evc,pfc,tspkc] = collspikes(xv,tv,pv,sv,spk); 
    ep = [ep mean(abs(evc))];  
end

