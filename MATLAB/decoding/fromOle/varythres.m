function varythres  

dt = 0.15; 
for thr=1:0.01:1.1
    fname = strcat('t',num2str(dt*1000));
    fname = strcat(fname,'thr');              
    fname = strcat(fname,num2str(100*thr));
    fname
    [tv,pv,ev,sv] = reconstruct4(dt,thr) ; 
    save(fname,'tv','pv','ev','sv');
end

