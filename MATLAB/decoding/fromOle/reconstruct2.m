function [xv,tv,pv,ev,sv] =reconstruct2(dt,div,p,fall,nall,occ)
%
% function [xv,tv,pv,ev,sv] =reconstruct2(dt,div,fall,nall,occ)
%
% dt  :  time window
% div :  # of phase divisions of the theta cycle
% pha :  0 = do not use phase advance (shuffled place fields) 
%        1 = use phase advance 
% cont:  0 = 1st half place fields, 2nd half reconstruct. 
%        1 = update place cells behind rat

bayesian = 1;
clf

%------------------------------------------------------------
% Define constants
%------------------------------------------------------------
bins = 150; 
lgd  = 205;
dts  = 0.05; 
cells = div*38; 
dev = 50; 

tmin = 1e4;               
tmax = 1.1e4;


%------------------------------------------------------------
% Allocate arrays         
%------------------------------------------------------------

xv = zeros(1,ceil((tmax-tmin)/dt));
tv = zeros(1,ceil((tmax-tmin)/dt));
pv = zeros(1,ceil((tmax-tmin)/dt));
ev = zeros(1,ceil((tmax-tmin)/dt));
sv = zeros(cells,ceil((tmax-tmin)/dt));


ploc = zeros(bins,1);




maxind = 0;
minind = 1e4;
for i=1:size(fall,2) 
     ind = find(fall(:,i) >= tmin &  fall(:,i) <= tmax);
     if min(ind) < minind
         minind = min(ind);
     end  
     if max(ind) > maxind
         maxind = max(ind);
     end  
end

fall = fall(minind:maxind,:);

ind = find(p(:,1) >=tmin &  p(:,1) <= tmax);
p = p(ind,:);


%------------------------------------------------------------
% Build place fields and reconstruct        
%------------------------------------------------------------


k = 0;
l = 0; 
ploc = occ./sum(occ); 
firstrun = 1; 
for t=(tmin+(tmax-tmin)/2):dt:tmax-2*dt          

    nt = getspikes(fall,t,dt);
    xt  = getbin(p,t+dt/2); 
    if firstrun
        firstrun = 0; 
        xold = xt;
    end
 
    if bayesian 
        for i=1:bins
           fx = (nall(i,:)/occ(i))';
           Pv = exp(-(geterror(i,xold)^2)/(2*(dev^2)));
           P(i) =  Pv.*prod(fx.^nt)*exp(-dt*sum(fx));
        end
    else
        for i=1:bins
           fx = (nall(i,:)/occ(i))';
           P(i) = nt'*fx;
        end
    end
 
    [dummy,xp] = max(P);
    if dummy ~= 0 
        err = geterror(xp,xt); 
        k = k + 1; 
        xold =  xp; 
        xv(k) = t; 
        tv(k) = xt;
        pv(k) = xp;
        ev(k) = err;
        sv(1:cells,k) = nt;        
   %     plot(xt,xp,'.'); 
        hold on
        if mod(k,10) == 0 
            fprintf(1,'Done %d Time=%d  Error=%d   Mean err=%d \n',(t-tmin)/(tmax-tmin),t,abs(err),mean(abs(ev(1:k)))); 
     %       drawnow;
     %       set(gca,'XLim',[0 151]);
     %       set(gca,'YLim',[0 151]);
        end

    end
end 
    
xv = xv(1:k);
pv = pv(1:k);
tv = tv(1:k);
ev = ev(1:k);
sv = sv(:,1:k); 
