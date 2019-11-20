function nallburst
% Contruct place fields (nall) and spike falls (fall)  appling the burst filter 

load pos
load occu
load ttheta

for dtb=0:0          
    fname = strcat('fallf',num2str(1000*dtb))  
    fall = cellfiring(dtb);                      
    save(fname,'fall'); 

    fname = strcat('nallf',num2str(1000*dtb))  
    nall = buildplacefield(p,ttheta,occu,1,dtb);
    save(fname,'nall'); 

end
