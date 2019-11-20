clear all
cd /home/ojensen/skaggs/reconstruct/jns

load pos
load occu
bayesian = 0;
phase    = 1; 
 


div = 1; 
for bayesian=1:1:1     
  for dt = 0.050:0.050:1.5               
    clear fallp fallr nallp
    fname = strcat('nall',num2str(div))  ;
    load(fname);
    fname = strcat('fall',num2str(div)) ;
    load(fname);


    fname = makename(dt,div,phase,bayesian)
    [xv,tv,pv,ev,sv] =reconstruct(dt,div,p,fallp,nallp,occu,bayesian);
    save(fname,'xv','tv','pv','ev','sv'); 
  end
end 

dt = 0.150;    
for bayesian=1:1:1     
  for div = 2:9                           
    clear fallp nallp
    fname = strcat('nall',num2str(div))  ;
    load(fname);
    fname = strcat('fall',num2str(div)) ;
    load(fname);


    fname = makename(dt,div,phase,bayesian)
    [xv,tv,pv,ev,sv] =reconstruct(dt,div,p,fallp,nallp,occu,bayesian);
    save(fname,'xv','tv','pv','ev','sv'); 
  end
end 
