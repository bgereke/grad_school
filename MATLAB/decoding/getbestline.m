function [best] = getbestline(Pxn)

%find longest contiguous path (time bins where some probability info exists)

z = reshape(Pxn,numel(Pxn),1);
cut = prctile(z,99);
%Pxn(Pxn<=cut) = 0;

%find peak prob space bin for each time bin
for i = 1:size(Pxn,2)
   
   [m,ind] = max(Pxn(:,i)); 
   if m>0
       maxinds(i,1) = ind;
       com(i,1) = round((Pxn(:,i)'*(1:1:length(Pxn(:,i)))')/sum(Pxn(:,i))); %compute com also
   else
       maxinds(i,1) = nan;
       com(i,1) = nan;
   end
end

%find longest contiguous line
ind = ~isnan(maxinds);
dif = diff([0 ind' 0]);
up = find(dif == 1);
dn = find(dif == -1);
[~, ind] = max(dn - up);
k = [up(ind) dn(ind)-1];%the index of the start and stop of longest line

if ~isempty(k)
best(:,1) = k(1,1):k(1,2);
best(:,2) = maxinds(k(1,1):k(1,2),1);
best(:,3) = com(k(1,1):k(1,2),1);
else
    best = [];
end

