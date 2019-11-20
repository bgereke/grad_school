function v = makevelo(p)
% function v = makevelo(p)
%
% p = [time x y angle]
% v = [time velocity]
% 
% Generate a vector of velocities from p 

v = zeros(length(p),1); 
l = zeros(length(p),1); 

for i=1:length(p)
    l(i) = getdist([p(i,2) p(i,3)]);
   if mod(i,1000) == 0
       fprintf(1,'%d ',floor(100*i/length(l)));
   end
end

for i=1:length(l)-1; 
   v(i) = (l(i+1)-l(i))/(p(i+1)-p(i));
   if mod(i,1000) == 0
       fprintf(1,'%d ',floor(100*i/length(l)));
   end
   if v(i) < -1000
      v(i) = v(i-1);
   end
end

v = [p(:,1) v];
