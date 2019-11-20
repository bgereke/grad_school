function p = locate(x,v)
% function p = locate(x,v)
%
% Find the indices of the item in x being less or equal 
% than v. Elements in x must increase in order. 

x = x(:,1);  
p2 = length(x);
p1 = 1; 
while p1 ~= p2
    pm = p1 + ceil((p2 - p1)/2);
    if v <  x(pm) 
        p2 = pm-1;                 
    else
        p1 = pm;       
    end
end
p = p1;

