function plotgamma(c1,c2,c3,c4,c5,c6,tmin,tmax); 



load c1_1
load c3_3
load c4_1
load c4_3
load tgamma 
load ttheta

subplot(6,1,1)
s = ones(size(c1));
stem(c1,s,'.')
set(gca,'XLim',[tmin tmax]); 
ylabel('t_{theta}')

subplot(6,1,2)
s = ones(size(c2));
stem(c2,s,'.')
set(gca,'XLim',[tmin tmax]); 
ylabel('t_{gamma}')

subplot(6,1,3)
s = ones(size(c3));
stem(c3,s,'.')
set(gca,'XLim',[tmin tmax]); 
ylabel('CA3, pyram')

subplot(6,1,4)
s = ones(size(c4));
stem(c4,s,'.')
set(gca,'XLim',[tmin tmax]); 
ylabel('CA3, inhi')

subplot(6,1,5)
s = ones(size(c5));
stem(c5,s,'.')
set(gca,'XLim',[tmin tmax]); 
ylabel('DG, pyram')

subplot(6,1,6)
s = ones(size(c6));
stem(c6,s,'.')
set(gca,'XLim',[tmin tmax]); 
ylabel('DG, inhi')
