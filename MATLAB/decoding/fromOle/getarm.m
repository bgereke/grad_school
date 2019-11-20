function a = getarm(l)
% function a = getarm(l)
% x = [x1 x2] 
%
% Return the arm (1,2 or 3) of the maze accroding to location x

y1 =  0.7179 * l(1) + 69;
y2 = -0.4528 * l(1) + 209;
y3 = -9 * l(1) + 1254;
% plot(x,[y1' y2' y3']); 

a = 1; 
if l(2) >= y2 & l(2) >= y1   
    a = 3;
end

if l(2) <= y3 & l(2) <= y2   
    a = 1;
end

if l(2) >= y3 & l(2) <= y1   
    a = 2;
end
