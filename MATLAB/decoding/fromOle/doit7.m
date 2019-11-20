function doit7

for dt=1:-0.05:0.05
    fname = strcat('pre',num2str(dt*100));
    fname
    [tv,pv,ev,sv] = reconstruct2(dt) ; 
    [tv2,pv2,ev2,sv2] = reconstruct3(dt)  ;
    save(fname,'tv','pv','ev','sv','tv2','pv2','ev2','sv2');
end


