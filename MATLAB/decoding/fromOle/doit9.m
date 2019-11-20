dt = 0.150; 
for div=1:9    
    fname = strcat('cd',num2str(div));
    fname = strcat(fname,'t');
    fname = strcat(fname,num2str(1000*dt));
fname
    [xv,tv,pv,ev,sv,nall,occ] =reconcont(dt,div);
    save(fname,'xv','tv','pv','ev','sv','nall','occ'); 
end


div = 1; 
for dt=0.050:0.050:1.500  
    fname = strcat('cd',num2str(div));
    fname = strcat(fname,'t');
    fname = strcat(fname,num2str(1000*dt));
fname
    [xv,tv,pv,ev,sv,nall,occ] =reconcont(dt,div);
    save(fname,'xv','tv','pv','ev','sv','nall','occ'); 

end
