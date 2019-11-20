clear all
cd /home/ojensen/skaggs/reconstruct

load pos
load occu
phase = 1; 
spk = 5;   

for div = 6:9                          
    clear fallp nallp
    fname = strcat('nall',num2str(div))  
    load(fname);
    fname = strcat('fall',num2str(div))  
    load(fname);

    fname = makename4(div,spk,phase)
    [xv,tv,pvb,pvt,evb,evt,sv] =reconstruct3(spk,div,p,fallp,nallp,occu);
    save(fname,'xv','tv','pvb','pvt','evb','evt','sv'); 

end



