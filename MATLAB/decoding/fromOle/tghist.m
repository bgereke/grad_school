function tghist(x,y);


Tmin = 0.100  
Tmax = 0.167

subplot(2,1,1); 
x = 0.0100*(rand(length(x),1)-0.5)+ x; 
z = x; 
j = 0; 
for i=1:length(x)
    if x(i) > Tmin & x(i) < Tmax   
        j = j + 1; 
        z(j) = x(i) ; 
    end 
    if mod(i,1000) == 0
        fprintf('%d\n',i/length(x)); 
    end
end
z  = 1000*z(1:j);  
hist(z,50) 
xlabel('T_{theta} (msec)') 

set(gca,'XLim',[1000*Tmin 1000*Tmax]); 

x =y ; 
subplot(2,1,2); 
Tmin = 0.0125; 
Tmax = 0.0333;  


x = 0.0100*(rand(length(x),1)-0.5)+ x; 
z = x; 
j = 0; 
for i=1:length(x)
    if x(i) > Tmin & x(i) < Tmax   
        j = j + 1; 
        z(j) = x(i) ; 
    end 
    if mod(i,1000) == 0
        fprintf('%d\n',i/length(x)); 
    end
end
z  = 1000*z(1:j);  
hist(z,50) 
xlabel('T_{gamma} (msec)') 

set(gca,'XLim',[1000*Tmin 1000*Tmax]); 
