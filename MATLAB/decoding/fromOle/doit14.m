function [eb,et] =doit14

div = 1;           
phase = 1; 

eb = [];
et = [];
spk = 5;
for div = 1:6                           
spk
    bayesian= 0; 
    fname = makename4(div,spk,phase)
    load(fname)
    et(div) = mean(abs(evt));
    eb(div) = mean(abs(evb));
    clear evt evb pvb pvt sv xv tv 
end


