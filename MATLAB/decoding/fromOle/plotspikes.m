function plotspikes(tv,pv,sv,peaks)
close
plot(tv(1,:),tv(2,:),'y-')
hold on
fv = sum(sv(2:size(sv,1),:)); 
for i=1:size(tv,2)
    if fv(i) > 3       
        plot(pv(1,i),pv(2,i),'gx')
    end 
end


hold on
for j=1:length(peaks)
    for i=1:size(sv,2)
        if sv(j+1,i) > 0
           plot(sv(1,i),peaks(j),'k.')
         %  fcircle(sv(1,i),peaks(j),0.1*sv(j+1,i),'k');
        end
    end
end

