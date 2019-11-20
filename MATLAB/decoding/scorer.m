function [score] = scorer(seq)

%make sure seq is already in temporal order

t = seq(:,1);
c = seq(:,2);

dt = diffmatrix(t);
dc = diffmatrix(c);

if length(find(dt(1,:)<0)) > 0 
    disp('error : not sorted');
    return;
end

dtc = dt.*dc;

dtc(dtc<0) = -1;
dtc(dtc>0) = 1;
% t
% dt
% c
% dc
% dtc
score = sum(sum(dtc));