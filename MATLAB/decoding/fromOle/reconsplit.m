clear all
cd /home/ojensen/skaggs/reconstruct/jns

load pos
load occu
bayesian = 0;
phase    = 1; 
 
div = 5; 
dt = 0.150; 

for i=1:div
    fname = strcat('nall',num2str(div))  ;
    load(fname);
    fname = strcat('fall',num2str(div)) ;
    load(fname);
    [fall,nall] = subdiv(fallp,nallp,div,i); 
    fname = makename2(dt,div,phase,i,bayesian)
    [xv,tv,pv,ev,sv] =reconstruct(dt,1,p,fall,nall,occu,bayesian);

    save(fname,'xv','tv','pv','ev','sv'); 
end  
