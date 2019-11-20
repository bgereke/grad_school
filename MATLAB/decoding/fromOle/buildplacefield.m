function x = buildplacefield(p,theta,div,tmin,tmax)
% function x = buildplacefield(p,theta,div,tmin,tmax)
%
% dtb = threshold for the burst filter (0 for no filter).
%
% Build the place fields for the cells defined by tetrode and cn

fprintf(1,'start: function buildplacefield, tmin=%d tmax=%d\n',tmin,tmax); 

tetrode = [1 1 1 1 2 2 3 3 3 3 3 3 3 4 4 5 5 6 6 6 6 6 7 7 7 7 7 8 8 8 8 8 8  11 11 11 11 11]; 
cn      = [1 3 6 7 2 3 1 2 4 5 6 7 8 1 5 1 2 3 4 5 6 7 1 2 6 7 8 1 2 3 4 5 10 1  2  4  6  9]; 

occ = occupancy(p,tmin,tmax); 

x = [];
fprintf(1,'Total # cells = %d \n',length(tetrode)); 
fprintf(1,'Cell '); 
for i=1:length(tetrode)
   fprintf(1,'%d ',i); 
%   [dummy,s] = burstfilter(loadcell(tetrode(i),cn(i)),dtb);
   s = loadcell(tetrode(i),cn(i));
   s =  s(find(s >tmin &  s < tmax));

   x = [x placefield(s,p,theta,occ,div)];
end


fprintf(1,'\n');                                   
fprintf(1,'stop:  function buildplacefield\n'); 
