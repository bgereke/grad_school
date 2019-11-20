function w = burstfilter(x,Arel);

w = 0*x; 

A = 0;
j = 0;
for i=2:length(x);
    A = calcium(A,x(i) - x(i-1));  
    if A > Arel  
        j = j + 1;
        w(j) = x(i); 
%        A = 0; 
    end
end
w = w(1:j);  

