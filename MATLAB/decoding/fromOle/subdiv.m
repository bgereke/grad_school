function [fallo,nallo] = subdiv(fall,nall,div,c)


fallo =[];
nallo = [];

for i=1:38   
    j = (i-1)*div + c; 
    fallo = [fallo fall(:,j)]; 
    nallo = [nallo nall(:,j)]; 
end
