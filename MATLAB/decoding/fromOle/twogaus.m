function twogaus

subplot(2,1,1) 
x= 1:0.1:65;   
y1 = 10*gausfunc(x,30,5); 
y2 = 10*gausfunc(x,43,5);        
plot(x,y1,'y',x,y2,'y',x,y1+y2,'k')
set(gca,'YLim',[0 1.8]); 
title('place field'); 

subplot(2,1,2) 
y1 = 10*gausfunc(x,36.5,5); 
y2 = 10*gausfunc(x,36.5,5);        
plot(x,y1,'y',x,y2,'y',x,y1+y2,'k')
set(gca,'YLim',[0 1.8]); 
title('place field shiftet 7.5 cm ahead'); 
