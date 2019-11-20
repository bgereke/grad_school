function plotres(z,tit,xmin,xmax);
clf        
xf = xmin:xmax;
subplot(4,1,1)
plot(z(:,1),z(:,2),'c.')
hold on
plot(z(:,1),z(:,2)+1,'c.')
% c1-1
% yf = -0.0491*(xf-115)+1.517;
% c6-3  
yf = -0.0451*(xf-195)+1.3020; 
% plot(xf,yf,'k-');
hold off
set(gca,'XLim',[xmin xmax]);
set(gca,'YLim',[0 2]);
ylabel('phase'); 
title(tit); 
  
fprintf(1,'velocity\n'); 
subplot(4,1,2)
p = polyfit(z(:,1),z(:,3),1)  
yf = polyval(p,130)  
yf = polyval(p,xf);
plot(z(:,1),z(:,3),'c.',xf,yf,'k-')

set(gca,'XLim',[xmin xmax]);
ylabel('velocity (cm/sec)'); 



fprintf(1,'gamma\n'); 
subplot(4,1,3)
p = polyfit(z(:,1),z(:,4),1)  
yf = polyval(p,130)  
yf = polyval(p,xf);
plot(z(:,1),z(:,4),'c.',xf,yf,'k-')
set(gca,'XLim',[xmin xmax]);
ylabel('# gamma per theta'); 

fprintf(1,'Ttheta\n'); 
subplot(4,1,4)
q = rand(length(z),1); 
z(:,5) = 0.005*(q-0.5) + z(:,5); 
p = polyfit(z(:,1),z(:,5),1)  
yf = polyval(p,130) 
yf = polyval(p,xf);
plot(z(:,1),1000*z(:,5),'c.',xf,1000*yf,'k-')
set(gca,'XLim',[xmin xmax]);
xlabel('position (cm)'); 
ylabel('T_{theta} (msec)'); 

orient tall
