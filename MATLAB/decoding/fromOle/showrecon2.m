function [lev1,lev2] = showrecon(dt)

load pos
load nallf 
load fall 
load nall4f
load fall4
load location 

showplot = 0;
ev1 = [];
ev2 = [];
bins = 150;
lgd = 205; 
sig = 20.0; 

perr1 = zeros(bins,2);
perr2 = zeros(bins,2);
lev1 = zeros(bins,1);
lev2 = zeros(bins,1);


tmin = p(1,1);
tmax = p(length(p),1) - 2*dt;

i = 0; 
if showplot
    clf 
end
pos  = getpos(p,tmin); x  = pos(1); y  = pos(2); l  = getdist([x y]);
xp1 = round(l*bins/lgd);
xp2 = xp1; 
for t=tmin:dt:tmax          
    [xp1, xt1, err1] = baysian(t,p,nallf,fall,loc,dt,sig,xp1);
    [xp2, xt2, err2] = baysian(t,p,nall4f,fall4,loc,dt,sig,xp2);
    if (xp1 ~= -1)
        ev1 = [ev1 err1];
        xt1r = floor(xt1);
        if xt1r < 1 
            xt1r = 1;
        end 
        if xt1r > bins
            xt1r = bins;
        end 
        perr1(xt1r,1) = perr1(xt1r,1) + err1;
        perr1(xt1r,2) = perr1(xt1r,2) + 1;
        for j=1:bins                 
            if perr1(j,2) ~= 0
                lev1(j) = perr1(j,1)./perr1(j,2);
            end
        end
    end
    if (xp1 ~= -1) & showplot
        subplot(2,1,1);
        hold on 
        plot(xt1,xp1,'r.'); 
    end

    if (xp2 ~= -1)
        ev2 = [ev2 err2];
        xt2r = floor(xt2);
        if xt2r < 1
            xt2r = 1;
        end
        if xt2r > bins 
            xt2r = bins;
        end
        perr2(xt2r,1) = perr2(xt2r,1) + err2;
        perr2(xt2r,2) = perr2(xt2r,2) + 1;
        for j=1:bins
            if perr2(j,2) ~= 0
                lev2(j) = perr2(j,1)./perr2(j,2);
            end  
        end
    end

    if (xp2 ~= -1) & showplot
        hold on 
        subplot(2,1,1);
        plot(xt2,xp2,'b.'); 
    end
    if mod(i,10) == 0 & showplot
        subplot(2,1,2);
        plot(1:bins,lev1,'r',1:bins,lev2,'b');
        drawnow;
        fprintf(1,'No phase %d With phase %d Improvement %d  \n',mean(ev1),mean(ev2),100*(mean(ev1)-mean(ev2))/mean(ev1))    
    end
    i = i + 1  ;
%        fprintf(1,'No phase %d With phase %d Improvement %d  \n',mean(ev1),mean(ev2),100*(mean(ev1)-mean(ev2))/mean(ev1))    
end 
e1 = mean(ev1);
e2 = mean(ev2);
