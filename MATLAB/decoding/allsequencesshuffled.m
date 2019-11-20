function [all] = allsequencesshuffled(spikes)
%this version returns all possible sequences after reassigning each cell a
%new set of COMs (from the existing ones)

ind = [];
for i = 1:size(spikes(:,1),1)
    if size(spikes{i,1},1) > 0
        ind = [ind ; i];
    end
end
r = randperm(size(ind,1));
newind = ind(r);

newspikes = cell (size(spikes));
newspikes(ind,1) = spikes(ind,1);
newspikes(ind,2) = spikes(newind,2);
spikes = newspikes; clear newpsikes;

numcenters = [];
for i = 1:size(spikes(:,1),1)
    if ~isempty(spikes{i,1})
        numcenters  = [numcenters ; [size(spikes{i,2},2) i]];    
    end
end

poss = prod(numcenters(:,1));
pm = 1;

for i = 1:size(numcenters(:,1),1)
   choices = numcenters(i,1); %number of choices
   tpm = pm;
   for j = 1:choices-1
       pm = [pm,tpm];   
   end
       numcol = size(pm,2);
       blocksize = numcol / choices;
    for j = 1:choices
       for k = 1:blocksize
            pm(i,(j-1)*blocksize+k) = j;
       end
    end    
end

%create all possible combinations using the pm (possibilty matrix)
%the number moving down each column of pm tells you which COM to use for
%that cell
all = [];
numcells = size(pm,1);
for i = 1:poss
    temp2 = [];
    temp = [];
    for j = 1:numcells
        temp = spikes{numcenters(j,2),1};
        temp(:,2) = spikes{numcenters(j,2),2}(1,pm(j,i));
        temp2 = [temp2;temp];
    end
    [~,inds] = sort(temp2(:,1));
    all (:,1,i) = temp2(inds,1);
    all (:,2,i) = temp2(inds,2);
end    

