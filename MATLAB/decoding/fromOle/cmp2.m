function [dst1,dst2] = cmp2(div1,div2)

fname = strcat('bayp_d',num2str(div1),'dt150r');
load(fname);
[xvc,pvc,tvc,evc,pfc,tspkc] = collspikes(xv,tv,pv,sv,10);
[evr,xvr,tvr,pvr]=removefoodstands(evc,xvc,tvc,pvc); 

fname = strcat('bayp_d',num2str(div1),'dt150p');
load(fname);
[xvc,pvc,tvc,evc,pfc,tspkc] = collspikes(xv,tv,pv,sv,10);
[evp,xvp,tvp,pvp]=removefoodstands(evc,xvc,tvc,pvc); 

m = 0; 
for j=1:length(tvr)
    for k=1:length(tvp)
        if tvr(j) == tvp(k)
            m = m + 1;
 %           dst1(m) = abs(abs(evr(j)) - abs(evp(k)));
        end
    end
end

fname = strcat('bayp_d',num2str(div2),'dt150r');
load(fname);
[xvc,pvc,tvc,evc,pfc,tspkc] = collspikes(xv,tv,pv,sv,10);
[evr,xvr,tvr,pvr]=removefoodstands(evc,xvc,tvc,pvc); 

fname = strcat('bayp_d',num2str(div2),'dt150p');
load(fname);
[xvc,pvc,tvc,evc,pfc,tspkc] = collspikes(xv,tv,pv,sv,10);
[evp,xvp,tvp,pvp]=removefoodstands(evc,xvc,tvc,pvc); 

m = 0; 
for j=1:length(tvr)
    for k=1:length(tvp)
        if tvr(j) == tvp(k)
            m = m + 1;
 %           dst2(m) = abs(abs(evr(j)) - abs(evp(k)));
        end
    end
end
mean(dst1)
mean(dst2)

[H,Ptt] = ttest2(dst1,dst2,0.05,-1)


if length(dst1) < length(dst2)
    [Prk, H] = ranksum(dst1,dst2,0.05)
else
    [Prk, H] = ranksum(dst2,dst1,0.05)
end
