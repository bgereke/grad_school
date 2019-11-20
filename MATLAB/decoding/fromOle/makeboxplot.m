function [pt,mt] =makeboxplot

div = 1;           
load pos
load occu
phase = 1; 

pt = [];
pb = [];
for spk = 1:10                          
spk
    bayesian= 0; 
    fname = makename3(bayesian,div,spk,phase)
    load(fname)
    pt = [pt [spk prctile2(diff(xv),50) prctile2(diff(xv),25) prctile2(diff(xv),75) prctile2(diff(xv),5) prctile2(diff(xv),95)]'];
    mt(spk) = mean(diff(xv));
end


errorbar(1:10,pt(3,:),pt(1,:),pt(5,:),'+')
