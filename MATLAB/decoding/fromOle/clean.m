function pc = clean(p)
j=0; 
pc = p; 
for i=1:length(p)
    if p(i,1) > 9092 & p(i,1) < 12738  
        j = j + 1; 
        pc(j,:) = p(i,:) ; 
    end 
    if mod(i,1000) == 0
        fprintf('%d\n',i/length(p)); 
    end
end
pc = pc(1:j-1,:); 

