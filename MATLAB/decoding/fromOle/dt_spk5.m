function w = dt_spk5            

w =[];
for dt=0.050:0.050:1.0  
    fname = strcat('t',num2str(dt*1000));
    fname
    load(fname);
    [pvc,tvc,evc,pfc,tspkc] = collspikes(tv,pv,sv,0); 

    w = [w [dt mean(abs(evc))]']   

end

