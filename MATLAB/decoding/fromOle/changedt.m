function doit             

thr = 0; 
for dt=0.050:0.050:1.0  
    fname = strcat('t',num2str(dt*1000));
    fname
    [tv,pv,ev,sv] = reconstruct4(dt,thr) ; 
    save(fname,'tv','pv','ev','sv');
end

