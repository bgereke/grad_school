function w = placecellphases;  

tetrode = [1 1 1 1 2 2 3 3 3 3 3 3 3 4 4 5 5 6 6 6 6 6 7 7 7 7 7 8 8 8 8 8 8  11 11 11 11 11]; 
cn      = [1 3 6 7 2 3 1 2 4 5 6 7 8 1 5 1 2 3 4 5 6 7 1 2 6 7 8 1 2 3 4 5 10 1  2  4  6  9]; 

load ttheta
clear w;
for i=1:length(tetrode)
    fprintf(1,'%d \n',i); 
    x = loadcell(tetrode(i),cn(i));       
    y = zeros(size(x)); 
    for j=1:length(x)
        y(j) = getphase(ttheta,x(j)); 
    end
    [h,z] = hist(y,50);
    w(i,:) = h(26:50); 
end

