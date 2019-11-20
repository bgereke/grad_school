function [xv,tv,pvb,pvt,evb,evt,sv] =reconstruct3(minspk,div,p,fall,nall,occ)
%
% function [xv,tv,pvb,pvt,evb,evt,sv] =reconstruct3(minspk,div,fall,nall,occ)
%
% minspk: min spk's per time window
% div :  # of phase divisions of the theta cycle
% pha :  0 = do not use phase advance (shuffled place fields) 
%        1 = use phase advance 
% cont:  0 = 1st half place fields, 2nd half reconstruct. 
%        1 = update place cells behind rat


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
dt = 0.010;

%------------------------------------------------------------
% Allocate arrays         
%------------------------------------------------------------

xv = zeros(1,ceil((tmax-tmin)/0.050));
tv = zeros(1,ceil((tmax-tmin)/0.050));
pvt = zeros(1,ceil((tmax-tmin)/0.05));
pvb = zeros(1,ceil((tmax-tmin)/0.05));
evt = zeros(1,ceil((tmax-tmin)/0.050));
evb = zeros(1,ceil((tmax-tmin)/0.050));
sv = zeros(cells,ceil((tmax-tmin)/0.050));


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


tstart = tmin+(tmax-tmin)/2;
nt = zeros(cells,1);
xold  = getbin(p,tmin+(tmax-tmin)/2); 
told = tmin+(tmax-tmin)/2;
for t=(tmin+(tmax-tmin)/2):dt:tmax-2*dt          

    nt = getspikes(fall,t,dt) + nt;
    if sum(nt) >= minspk          
        % reconstruct  
        xt  = getbin(p,t+(t-tstart)/2); 

        %  Bayesian method 
        for i=1:bins
            fx = (nall(i,:)/occ(i))';
            Pv = exp(-(geterror(i,xold)^2)/(2*(dev^2)));
            P(i) =  Pv.*prod(fx.^nt)*exp(-(t-told)*sum(fx));
        end
        [dummy,xpb] = max(P);
        errb = geterror(xpb,xt); 
        %  Template method  
        for i=1:bins
            fx = (nall(i,:)/occ(i))';
            P(i) = nt'*fx;
        end 
        [dummy,xpt] = max(P);
        errt = geterror(xpt,xt); 
   

        k = k + 1; 
        xold =  xpb; xv(k) = t; tv(k) = xt; 
        told = t;
        pvb(k) = xpb; 
        pvt(k) = xpt; 
        evb(k) = errb; 
        evt(k) = errt; 

        sv(1:cells,k) = nt;        
%        subplot(2,1,1) 
%        plot(xt,xp,'.'); 
%        title(makename3(bayesian,div,minspk,1))
%        hold on
%        subplot(2,1,2) 
        if k > 110
%             hist(diff(xv(1:k)),100)
        else
             if k > 3
%                 hist(diff(xv(1:k)))
             end
        end
        if mod(k,10) == 0 
            fprintf(1,'Done %d Time=%d  Mean template err=%d Mean Bayesian err=%d dt=%d \n',(t-tmin)/(tmax-tmin),t,mean(abs(evt(1:k))),mean(abs(evb(1:k))),t-tstart); 
%            drawnow;
%            set(gca,'XLim',[0 151]);
%            set(gca,'YLim',[0 151]);
        end



        nt = zeros(cells,1);
        tstart = t;
    end

end


xv = xv(1:k);
pvt = pvt(1:k);
pvb = pvb(1:k);
tv = tv(1:k);
evt = evt(1:k);
evb = evb(1:k);
sv = sv(:,1:k); 
