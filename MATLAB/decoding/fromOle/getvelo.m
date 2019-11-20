function x = getvelo(v,t)
% function x = getvelo(v,t)
% v = [time velocity]                  
% x = velocity        
% Find the velocity at time t                                      
%
% Ole Jensen, May 26, 1998

if size(v,2) > size(v,1)
    v = v';
end


j = locate(v,t); 
x = v(j,2);  

