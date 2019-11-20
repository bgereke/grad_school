function  [dv,err] =  divcoll

dv = [];
err= [];
for d=10:11  
    fname = 't150'
    fname = strcat(fname,'d');           
    fname = strcat(fname,num2str(d));
    fname  
    load(fname)  

    [pvc,tvc,ev,pfc,tspkc] = collspikes(tv,pv,sv,0); 
    dv = [dv d];
    err = [err mean(abs(ev))]  

end
