function x = loadcell(t,n)
fname = strcat('c',num2str(t));
fname = strcat(fname,'_'); 
fname = strcat(fname,num2str(n));

fid = fopen(fname,'r');
x = fscanf(fid,'%e');
fclose(fid); 
