function w = doit(nspk); 
 % clear all
cd /home/ojensen/skaggs/reconstruct/jns

load pos
load occu
bayesian = 0;
phase    = 1; 
 
div = 5; 
dt = 0.150; 
w = [];
for i=1:div
    fname = makename2(dt,div,phase,i,bayesian)
    load(fname)
    [xvc,pvc,tvc,evc,pfc,tspkc] = collspikes(xv,tv,pv,sv,nspk);

    w = [w [i mean(evc) length(xvc) ]'];
end  
