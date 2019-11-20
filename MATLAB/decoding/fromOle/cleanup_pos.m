function z = cleanup_pos(x)
% function z = cleanup_pos(x)
% Removes unreaden positions                             

z = x; 
j=0; 
for i=1:length(x)
    if x(i,2) ~= 0
        j = j+1; 
        z(j,:) = x(i,:);
    end
end
z=z(1:j,:);

% z(:,1) = z(:,1) - z(1,1); 
