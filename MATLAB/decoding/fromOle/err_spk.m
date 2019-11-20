function w = err_spk

w = [];
for i=75:5:100
    fname = strcat('pre',num2str(i)); 
    load(fname) 
    [pvc,tvc,evc]    = collspikes(tv,pv,sv,5);
    [pvc2,tvc2,evc2] = collspikes(tv2,pv2,sv2,5);
    w = [w [i 205*mean(abs(evc))/150  205*mean(abs(evc2))/150  205*median(abs(evc))/150  205*median(abs(evc2))/150]'];
end 
