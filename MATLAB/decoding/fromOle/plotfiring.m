function plotfiring(f,w,p)
% function plotfiring(f,w,p)
%  
%  f: firing of cell  
%  w: theta phase
%  p: position
%
clf   
plot(p(:,1),p(:,2),'k.'); 
hold on 
y = zeros(length(f),3); 
for i=1:length(f)
    l = getpos(p,f(i)); 
    r = getphase(w,f(i))  ;
    r = mod(r+0.5,1); 
    if r >= 0 & r <= 0.33
       plot(l(1),l(2),'r.'); 
    end 
    if r > 0.33 & r < 0.67
       plot(l(1),l(2),'y.'); 
    end 
    if r >= 0.67 & r <= 1
       plot(l(1),l(2),'b.'); 
    end 
    y(i,1:3) = [l(1) l(2) r];
    if mod(i,100) == 0
        fprintf(1,'%d\n',i/length(f)); 
    end
end
p1_1 = y; 
% save p1_1.mat p1_1  

 
