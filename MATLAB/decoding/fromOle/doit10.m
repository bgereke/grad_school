dt = 0.150; 
for div=1:4:5  
    fname = strcat('cpd',num2str(div));
    fname = strcat(fname,'t');
    fname = strcat(fname,num2str(1000*dt));
fname
    [xv,tv,pv,ev,sv,nall,occ] =reconcont2(dt,div);
    save(fname,'xv','tv','pv','ev','sv','nall','occ'); 
end


