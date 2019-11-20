function w = errajumps(fname)

load(fname)

for i=0:15
    [xvc,pvc,tvc,evc,pfc,tspkc] = collspikes(xv,tv,pv,sv,i);
    err = mean(abs(evc));
    jmp = length(find(abs(evc > 10 )))/length(evc);    
    w = [w [i err jmp]']  
end


  
