function x = cellfiring(tmin,tmax)
% function x = cellfiring(tmin,tmax)
%
% x = [spike-times x cells]
% Returns the spike times for the cells defined by tetrode
% and cn. The rows of spike times for each cell are ordered 
% are sorted by increasing time. The number of rows in x are 
% determined by the cell with % the maximum # of spikes. 
% Shorter columns are padded with 'inf' 



tetrode = [1 1 1 1 2 2 3 3 3 3 3 3 3 4 4 5 5 6 6 6 6 6 7 7 7 7 7 8 8 8 8 8 8  11 11 11 11 11]; 
cn      = [1 3 6 7 2 3 1 2 4 5 6 7 8 1 5 1 2 3 4 5 6 7 1 2 6 7 8 1 2 3 4 5 10 1  2  4  6  9]; 

lgd = 0;
for i=1:length(tetrode)
   s = loadcell(tetrode(i),cn(i));
   s =   s(find(s >tmin &  s < tmax));

   if (length(s) > lgd) 
       lgd = length(s);
   end
end

x = Inf*ones(lgd,length(cn));

for i=1:length(tetrode)
%   [dummy,s] = burstfilter(loadcell(tetrode(i),cn(i)),dtb);
   s = loadcell(tetrode(i),cn(i));
   s =   s(find(s >tmin &  s < tmax));
   fprintf(1,'%d ',i); 
   x(1:length(s),i) = s;
end

