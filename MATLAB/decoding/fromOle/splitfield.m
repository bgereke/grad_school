function [wsp,wsr]  = splitfield(w,theta,div)
% function [wsp,wsr] = splitfield(w,theta,div)
%
% Divide spikes times into 'div' divisions according to theta phase.
% If theta phase for a given spike time is -1 (not defined) the spikes
% is removed. w and ws are padded with Inf for spike vectors shorther
% than max length.

% w     : [spike-time cells]  
% theta : [theta-peaks]
% div   : number of division   
% wsp   : [spike-time cells]  
% wsr   : [spike-time cells]  


wsp = Inf*ones(size(w,1),size(w,2)*div);
wsr = Inf*ones(size(w,1),size(w,2)*div);
posp = zeros(size(w,2)*div);
posr = zeros(size(w,2)*div);

lgd = size(w,1); 

for i=1:size(w,2)
    fprintf(1,'%d ',i);
    for k=1:lgd                   
        if w(k,i) ~= Inf
            phase = getphase(theta,w(k,i));
            if phase == 1
                phase = 0.999;
            end
            if phase ~= -1
                c  = 1+floor(phase*div);
                l = c + (i-1)*div; 
                posp(l) = posp(l) + 1;
                wsp(posp(l),l) = w(k,i); 

                c =  1+floor(div*rand(1,1));
                l = c + (i-1)*div; 
                posr(l) = posr(l) + 1;
                wsr(posr(l),l) = w(k,i); 


            end
        end
    end
end

