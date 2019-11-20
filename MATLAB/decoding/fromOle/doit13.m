function [divc,evv] = doit13;
dt = 0.150;
load pos
load occu
evv = [];

for div = 1:9      
    fname = makename(dt,div,1,0);             
    load(fname);
    [xvc,pvc,tvc,evc,pfc,tspkc] = collspikes(xv,tv,pv,sv,10);       
    mean(abs(evc))
    evv = [evv mean(abs(evc))]; 
end

divc = 1:9;
