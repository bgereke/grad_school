function varydiv 

dt = 0.15; 
for div=10:11
    fname = strcat('t',num2str(dt*1000));
    fname = strcat(fname,'d');              
    fname = strcat(fname,num2str(div));
    fname
    [tv,pv,ev,sv] = reconstruct3(dt,div) ; 
    save(fname,'tv','pv','ev','sv');
end

