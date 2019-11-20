function [w1,w2] = rm_foodstands(x1,x2)


s1 = 1; 
s2 = 75.7;
s3 = 150.0;
wd = 3; 
w1 = 0*x1;
w2 = 0*x2;
j = 0; 
for i=1:length(x1)
    if ((x1(i) > s1+wd) & (x1(i) < s2-wd))  | ((x1(i) > s2+wd) & (x1(i) < s3-wd))
        j = j + 1;
        w1(j) = x1(i); 
        w2(j) = x2(i); 
    end
end
w1 = w1(1:j); 
w2 = w2(1:j); 
