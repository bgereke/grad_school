function showcorr(x)
% x = [time theta ratio velocity]

subplot(2,2,1)
x(:,2) = 0.0050*(rand(length(x),1)-0.5)+ x(:,2); 
p = polyfit(x(:,4),1000*x(:,2),1)
x2 = min(x(:,4)) : max(x(:,4));
y2 = polyval(p,x2); 
plot(x(:,4),1000*x(:,2),'c.',x2,y2,'k-')
% plot(x2,y2,'-')
xlabel('velocity (cm/sec)')
ylabel('T_{theta} (msec)')
s = num2str(p(1)); 
s = strcat(s,' v + '); 
s = strcat(s,num2str(p(2))); 
title(s)

subplot(2,2,2)
p = polyfit(1000*x(:,2),x(:,3),1)
x2 = min(1000*x(:,2)) : max(1000*x(:,2));
y2 = polyval(p,x2); 
plot(1000*x(:,2),x(:,3),'c.',x2,y2,'k-')
xlabel('T_{theta} (msec)')
ylabel('T_{theta}/T_{gamma}')
s = num2str(p(1)); 
s = strcat(s,' T_{theta}+ '); 
s = strcat(s,num2str(p(2))); 
title(s)

subplot(2,2,3)
p = polyfit(1000*x(:,2),1000*x(:,2)./x(:,3),1)
x2 = min(1000*x(:,2)) : max(1000*x(:,2));
y2 = polyval(p,x2); 
plot(1000*x(:,2),1000*x(:,2)./x(:,3),'c.',x2,y2,'k-')
xlabel('T_{theta} (msec)')
ylabel('T_{gamma} (msec)')
s = num2str(p(1)); 
s = strcat(s,' T_{theta}+ '); 
s = strcat(s,num2str(p(2))); 
title(s)

subplot(2,2,4)
p = polyfit(x(:,4),1000*x(:,2)./x(:,3),1)
x2 = min(x(:,4)) : max(x(:,4));
y2 = polyval(p,x2); 
plot(x(:,4),1000*x(:,2)./x(:,3),'c.',x2,y2,'k-')
xlabel('velocity (cm/sec)')
ylabel('T_{gamma} (msec)')
s = num2str(p(1)); 
s = strcat(s,' v + '); 
s = strcat(s,num2str(p(2))); 
title(s)

orient tall 
