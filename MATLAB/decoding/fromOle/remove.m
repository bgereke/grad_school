function y = remove(x,p)
% Remove location p in vector x and replaces the last
% element with 'inf'

y =[x(1:p-1,:)' x(p+1:length(x),:)']';
y(length(x),:) = Inf*ones(1,size(x,2));        
