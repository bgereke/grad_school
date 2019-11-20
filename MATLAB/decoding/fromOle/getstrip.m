function y=getstrip(x,xmin,xmax)

vmin = -10.0; 
y = zeros(size(x)); 
j = 0;
for i=1:length(x)
    if x(i,1) > xmin & x(i,1) < xmax & x(i,3) > vmin & x(i,2) > 0 
        j = j + 1;
        y(j,:) = x(i,:);  
    end
end
y = y(1:j,:); 
