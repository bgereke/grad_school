function [fq] = greedyBP(T_f,R)

fq = zeros(1,size(T_f,2));

for i = 1:size(T_f,2)
    
    [~,fq(i)] = max(T_f(:,i));
%     range = idx-R:idx+R;
%     cs = cumsum(T_f(range,i))/max(cumsum(T_f(range,i)));
%     fq(i) = find(cs>0.5,1,'first'); %median
    
end