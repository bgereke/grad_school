clear all
cd /home/ojensen/skaggs/reconstruct/jns

load pos
load occu
phase    = 0; 
 


dt = 0.150; 
for bayesian=1:1:1     
  for div =2:3                           
  %  clear fallp fallr nallp
    fname = strcat('nall',num2str(div))  ;
    load(fname);
    fname = strcat('fall',num2str(div)) ;
    load(fname);

    fname = makename(dt,div,phase,bayesian)
    [xv,tv,pv,ev,sv] =reconstruct(dt,div,p,fallr,nallp,occu,bayesian);
    save(fname,'xv','tv','pv','ev','sv'); 
  end
end 
