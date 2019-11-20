function [f3,m3] = localfit_ad(t,x,tmax,error,n)

tmax = tmax/2;
fit1 = zeros(size(x));
m1 = nan(size(x));
t = reshape(t,length(t),1);  
x = reshape(x,length(x),1);

for i = 1:length(t)
    support = t(abs(t(i)-t)<=tmax);    
    targets = x(abs(t(i)-t)<=tmax);
    if isnan(x(i))
        fit1(i) = nan;
    elseif sum(~isnan(targets)) == 1
        fit1(i) = x(i);
%     elseif sum(~isnan(targets)) == 2
%         tar = targets(~isnan(targets));
%         sup = support(~isnan(targets));
%         fit1(i) = x(i);
%         m1(i) = (tar(2)-tar(1))/(sup(2)-sup(1));
    else
        tar = targets(~isnan(targets));
        sup = support(~isnan(targets));
        tidx = find(sup==t(i));
        sdiffs = abs(t(i)-sup);
        [~,addidx] = sort(sdiffs,'ascend');  
%         fit1(i) = x(i);
%         m1(i) = (tar(addidx(1))-tar(addidx(2)))/(sup(addidx(1))-sup(addidx(2)));
%         if addidx(1)<addidx(2)
%             m1(i) = -1*m1(i);
%         end
        goodfit = 1; j = 1; fails = 0;
        ss = []; tt = [];
        while goodfit
            if j <= length(addidx) 
                ss = [ss; sup(addidx(j))];
                tt = [tt; tar(addidx(j))];
                j = j+1;           
            else
                goodfit = 0;
                continue
            end
            %solve regression using QR decomposition
            [Q,R] = qr([ones(size(ss)),ss],0);
            b = R\(Q'*tt);
            if sum(abs(tt-b(2)*ss-b(1)) > error) > 0
               if fails                   
                   goodfit = 0;
               end
               fails = 1;
               continue
            end
            fit1(i) = b(2)*t(i)+b(1);
            m1(i) = b(2); 
            fails = 0;
        end
        if fit1(i) == 0
            fit1(i) = tar(tidx);
            if tidx>1 && tidx < length(tar)
                [m, midx] = min([sup(tidx)-sup(tidx-1) sup(tidx+1)-sup(tidx)]);
                m1(i) = (tar(tidx-1+midx)-tar(tidx-2+midx))/m;
            elseif tidx > 1
                m1(i) = (tar(tidx)-tar(tidx-1))/(sup(tidx)-sup(tidx-1));
            else
                m1(i) = (tar(tidx+1)-tar(tidx))/(sup(tidx+1)-sup(tidx));
            end
        end
    end
end
m2 = nan(size(m1));
f2 = nan(size(fit1));
for i = 1:length(t)
    support = t(abs(t(i)-t)<=tmax);    
    mtargets = m1(abs(t(i)-t)<=tmax);
    mtar = mtargets(~isnan(mtargets));
    msup = support(~isnan(mtargets));
%     mtidx = find(msup==t(i));
    sdiffs = abs(t(i)-msup);
    [~,addidx] = sort(sdiffs,'ascend');
    
    if isnan(m1(i))
        m2(i) = nan;
    elseif length(addidx)<=n
        m2(i) = median(mtar(addidx));
    else
        m2(i) = median(mtar(addidx(1:n)));
    end
    
    ftargets = fit1(abs(t(i)-t)<=tmax);
    ftar = ftargets(~isnan(ftargets));
    fsup = support(~isnan(ftargets));
%     ftidx = find(fsup==t(i));
    sdiffs = abs(t(i)-fsup);
    [~,addidx] = sort(sdiffs,'ascend');
    
    if isnan(m1(i))
        f2(i) = nan;
    elseif length(addidx)<=n
        f2(i) = median(ftar(addidx));
    else
        f2(i) = median(ftar(addidx(1:n)));
    end
end

m3 = nan(size(m1));
f3 = nan(size(fit1));
for i = 1:length(t)
    support = t(abs(t(i)-t)<=tmax);    
    mtargets = m2(abs(t(i)-t)<=tmax);
    mtar = mtargets(~isnan(mtargets));
    msup = support(~isnan(mtargets));
%     tidx = find(sup==t(i));
    sdiffs = abs(t(i)-msup);
    [~,addidx] = sort(sdiffs,'ascend');
    
    if isnan(m1(i))
        m3(i) = nan;
    elseif length(addidx)<=n
        m3(i) = mean(mtar(addidx));
    else
        m3(i) = mean(mtar(addidx(1:n)));
    end
    
    support = t(abs(t(i)-t)<=tmax);    
    ftargets = f2(abs(t(i)-t)<=tmax);
    ftar = ftargets(~isnan(ftargets));
    fsup = support(~isnan(ftargets));
%     tidx = find(sup==t(i));
    sdiffs = abs(t(i)-fsup);
    [~,addidx] = sort(sdiffs,'ascend');
    
    if isnan(m1(i))
        f3(i) = nan;
    elseif length(addidx)<=n
        f3(i) = mean(ftar(addidx));
    else
        f3(i) = mean(ftar(addidx(1:n)));
    end
end
