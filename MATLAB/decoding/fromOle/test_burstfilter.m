function test_burstfilter(t1,t2,A)

load c1_1
x = c1_1;

w = burstfilter(x,A); 

subplot(2,1,1);
stem(x,ones(1,length(x)),'.') 
set(gca,'XLim',[t1 t2]);  
subplot(2,1,2);
stem(w,ones(1,length(w)),'.') 
set(gca,'XLim',[t1 t2]);  
