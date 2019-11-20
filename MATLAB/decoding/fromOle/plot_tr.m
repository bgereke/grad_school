function w = plot_tr(tv,pv,sv,min_spk)

x = zeros(1,length(tv)) ;  
y = zeros(1,length(tv));   
j = 0; 
ts = sum(sv(2:size(sv,1),:),1);
size(ts)  
for i=1:length(ts)
   if ts(i) >= min_spk 
       x(i) = 205*tv(2,i)/150;
       y(i) = 205*pv(2,i)/150;
   else
       j = j + 1;
       x(i) = -1;
       y(i) = -1;
   end
end
plot(x,y,'k.'); 
xlabel('True Position (cm)');
ylabel('Reconstructed Position (cm)');
set(gca,'XLim',[0 205]);
set(gca,'YLim',[0 205]);
axis square
w = [x' y'];
