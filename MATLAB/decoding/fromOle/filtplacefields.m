function x = filtplacefields(nall)

fprintf(1,'Start: filtplacefields\n'); 

x = nall;
fprintf(1,'Cell ');
for i=1:size(nall,2)
    fprintf(1,' %d',i);
    x(:,i) = smoothing(nall(:,i),2)';
end

fprintf(1,'Stop:  filtplacefields\n'); 
