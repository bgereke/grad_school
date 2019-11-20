function z = getvelo_ratio(v,r)
% function z  getvelo_ratio(v,r)
%
% v = [time velocity]                     
% r = [time theta ratio]                     
%
% z = [time theta ratio velocity]


w = zeros(length(r),1); 
for i=1:length(r)    
   w(i) = getvelo(v,r(i,1));

    if mod(i,100) == 0
        fprintf(1,'%d\n',i/length(r));
    end
end

z = [r w] ; 


