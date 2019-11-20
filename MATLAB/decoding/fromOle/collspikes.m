function [xvc,pvc,tvc,evc,pfc,tspkc] = collspikes(xv,tv,pv,sv,nspk)
% function [xvc,pvc,tvc,evc,pfc,tspkc] = collspikes(xv,tv,pv,sv,nspk)
%
% Cleans the true and estimated position vectors from time-windows
% having less than nspk. Calculates position error and/or spiking
% statistics.  
%
% tv    : [time pos],   
%         true location
% pv    : [time pos],    
%         predicted location
% sv    : [time spk1 spk2 ... spkn]
%         # of spikes for individual place cells    
% nspk  : time windows with less than nspk will be cleaned out    
% 
% 
% 
% tv    : [time pos],   
%         true location, cleaned 
% pv    : [time pos],    
%         predicted location,cleaned  
% ev    : [time error]
%         the errors calculated from tv and pv  
% pfc   : [time number_of_cells]
%         # of cells that fires 
% tspkc : [time number_of_spikes]               
%         total # of spikes               

bins = 150;
evc = zeros(size(xv));
xvc = zeros(size(xv));
pfc = zeros(1,size(tv,2));
pvc = zeros(size(tv));
tvc = zeros(size(pv)); 
tspkc = zeros(size(xv));             
j = 0;
for i=1:length(tv)

    xt = tv(i);
    xp = pv(i);
    spks = sum(sv(1:size(sv,1),i)); 
    tmp  = sv(1:size(sv,1),i);
    pf = size(find(tmp'>0),2);
    tspk = sum(sv(1:size(sv,1),i));

    if spks >= nspk
        j = j + 1;   
        evc(j) = geterror(xp,xt);
        xvc(j) = xv(i);
        tvc(j) = tv(i);
        pvc(j) = pv(i);
        pfc(j) = pf;             
        tspkc(j) = tspk;             
        
    end

end
j
evc = evc(1:j);
xvc = xvc(1:j);
tvc = tvc(1:j); 
pvc = pvc(1:j); 
pfc = pfc(1:j); 
tspkc = tspkc(1:j); 

