function z = extract(f,w,p,v,r)
% function extract(f,w,p,v,r)
%
% f = [time of spike]
% w = [time of theta peaks]
% p = [time xpos ypos angle]
% v = [time velocity]                     
% r = [time theta ratio]                     
%
% z = [d phase velo ratio]; 


a(1) = 143.3;
a(2) = 23.8;
a(3) = 264.0;
tol = 45; 
z = zeros(length(f),5); 
for i=1:length(f)    
   phase = getphase(w,f(i)); 
   x     = getpos(p,f(i)); 
   d     = getdist([x(1) x(2)]); 
   arm   = getarm([x(1) x(2)]); 
   if x(3) < a(arm)-tol | x(3) > a(arm)+tol
           phase = -1;
   end 
   velo = getvelo(v,f(i));
   ratio = getratio(r,f(i));
   z(i,1:5) = [d phase velo ratio(1) ratio(2)]; 
   if mod(i,100) == 0
      fprintf('%d\n',i/length(f)); 
   end
end



