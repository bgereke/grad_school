function [fields,spikes,ratesegments,xbins] = getfields(cellist,binsize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load stuff
scale = .275;
[x,y,t] = readVideoData('VT1.nvt',scale);
[x,y] = axesRotater(x,y);
x = x - min(x);
v = findVelLinear(x);
[xL] = linearizedirection(x,t,x,t,smooth(v,15));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = 0;
ff = 0;
fid=fopen(cellist);
while 1
       tline = fgetl(fid);
       if ~ischar(tline) 
           break
       end
       i = i+1;
       spikes{i,1} = loadSpikes2(tline);
       spikes{i,2} = getSpikePos(spikes{i,1},x,y,t);
       spikes{i,3} = linearizedirection(spikes{i,2},spikes{i,1},x,t,v);
       [spikes{i,4}, xbins] = ratemap_decode(xL,spikes{i,3},binsize);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get basic stuff
numcells = size(spikes,1);
numbins = length(xbins);
maxrate = 0;
for i = 1:numcells
    temp = max(spikes{i,4});
    if temp > maxrate
        maxrate = temp;
    end
end
minaccept = .01 * maxrate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find the idealized fields
cells = 0;
for i = 1:numcells 
    temp = spikes{i,4}; %temp rate map
    
    ind = find(temp>=minaccept);
    if length(ind) > 2
        jumps = find(diff(ind)>2); %when there is a break in the contiguity
        
        fieldind = [];
        if ~isempty(jumps)
            fieldind(1,1) = ind(1);
            for j = 1:length(jumps)
                fieldind(j,2) = ind(jumps(j));
                fieldind(j+1,1) = ind(jumps(j)+1);
            end
                fieldind(j+1,2) = ind(end);
        else
            fieldind(1,1) = ind(1);
            fieldind(1,2) = ind(end);
        end
        
        %copy the field segments
        segments = 0;
        ratesegments = cell(1,1);
        for j = 1:size(fieldind,1)
            if fieldind(j,2) - fieldind(j,1) > -1
                segments = segments + 1;
                ratesegments{segments,1} = temp (fieldind(j,1):fieldind(j,2));
                ratesegments{segments,2} = (fieldind(j,1):fieldind(j,2));                                
            end
        end
        
        %finally cpoy to cells
        %make sure at least one field has a field with 3 contiguous >2Hz
        usecell = 0;
        if ~isempty(segments)
            for j = 1:segments
                tdiff = diff(find(ratesegments{j,1} >= 2));
                ttdiff = find(tdiff == 1);
                tttdiff = diff(ttdiff);
                isthere = length(find(tttdiff == 1));
                if isthere > 0
                    usecell = 1;
                end
            end
        
            
            if usecell == 1
                cells = cells + 1;
                vv = 0;
                for j = 1:segments
                    %if sum(ratesegments{j,1}) > 15 %condition to include the rate segment
                    vv = vv+1;
                    fields{cells,1}(:,vv) = zeros(numbins,1);
                    fields{cells,1}(ratesegments{j,2},vv) = ratesegments{j,1};
                    %get COM for the segment "field"
                    com = dot(fields{cells,1}(:,vv),  1:1:numbins) / sum(fields{cells,1}(:,vv));
                    fields{cells,2}(1,vv) = com;
                    fields{cells,3} = spikes{i,3};
                    fields{cells,4} = spikes{i,4};
                    fields{cells,5} = spikes{i,1};
                    %end
                end
            end
        end
 
    end
end
