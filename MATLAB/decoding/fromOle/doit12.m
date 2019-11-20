
cd /home/ojensen/skaggs/reconstruct
clear all

load pos
load occu
phase = 1; 
spk = 10; 

for div = 2:9                          
    fname = strcat('nall',num2str(div))  
    load(fname);
    fname = strcat('fall',num2str(div))  
    load(fname);

    bayesian= 0; 
    fname = makename3(bayesian,div,spk,phase)
    [xv,tv,pv,ev,sv] =reconstruct3(spk,div,p,fallp,nallp,occu,bayesian);
    save(fname,'xv','tv','pv','ev','sv'); 

    bayesian= 1; 
    fname = makename3(bayesian,div,spk,phase)
    [xv,tv,pv,ev,sv] =reconstruct3(spk,div,p,fallp,nallp,occu,bayesian);
    save(fname,'xv','tv','pv','ev','sv'); 
end



