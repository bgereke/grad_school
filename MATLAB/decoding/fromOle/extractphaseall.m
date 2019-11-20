function extractphaseall;

load ttheta
load pos

tetrode = [1 1 1 1 2 2 3 3 3 3 3 3 3 4 4 5 5 6 6 6 6 6 7 7 7 7 7 8 8 8 8 8 8  11 11 11 11 11]; 
cn      = [1 3 6 7 2 3 1 2 4 5 6 7 8 1 5 1 2 3 4 5 6 7 1 2 6 7 8 1 2 3 4 5 10 1  2  4  6  9]; 

fprintf(1,'Total # cells = %d \n',length(tetrode)); 
fprintf(1,'Cell '); 
for i=1:length(tetrode)
   fprintf(1,'%d ',i); 
   s = loadcell(tetrode(i),cn(i));
   w = extractphase(s,ttheta,p); 
   fname = strcat('phase',num2str(i));
   save(fname,'w'); 
end


