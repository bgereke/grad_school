function layouts = layoutList

mfile = mfilename('fullpath');
path = fileparts(mfile);
fid = fopen([path filesep '..' filesep '..' filesep 'defaults' filesep 'defaultLayout.txt']);
txt = textscan(fid,'%s','delimiter','\t');
fclose(fid);
txt = txt{1};
txt = reshape(txt,37,length(txt)/37)';
fld = txt(1,2:4:end);
names = txt(3:end,1);

layouts = {};
groups{1} = {'variablesTable','experimentProperties','worldsMenu'};
groups{2} = {'objectProperties','worldSketch','worldDrawing'};
groups{3} = {'shapeProperties','textureSketch','textureDrawing'};
groupNum = zeros(1,length(fld));
for ndx = 1:length(fld)
    for g = 1:length(groups)
        for h = 1:length(groups{g})
            if strcmp(groups{g}{h},fld{ndx})
                groupNum(ndx) = g;
            end
        end
    end
end
for ndx = 1:length(names)
    str = names(ndx);
    for f = 1:length(fld)
        pos = txt(ndx+2,4*(f-1)+2:4*(f-1)+5);
        if ~strcmp(pos{1},'-')
            pos = cellfun(@(x)str2double(x),pos);
            str{end+1} = fld{f}; %#ok<AGROW>
            str{end+1} = pos; %#ok<AGROW>
        end
    end
    if strcmp(names{ndx},'Default')
        layouts{end+1} = addL([{'Experiment'} str(reshape([2*find(groupNum==1); 2*find(groupNum==1)+1],1,2*length(find(groupNum==1))))]); %#ok<AGROW>
        layouts{end+1} = addL([{'World'} str(reshape([2*find(groupNum==2); 2*find(groupNum==2)+1],1,2*length(find(groupNum==2))))]); %#ok<AGROW>
        layouts{end+1} = addL([{'Texture'} str(reshape([2*find(groupNum==3); 2*find(groupNum==3)+1],1,2*length(find(groupNum==3))))]); %#ok<AGROW>
    else
        layouts{end+1} = addL(str); %#ok<AGROW>
    end
end


function layout = addL(inp)

layout.name = inp{1};
for ndx = 2:2:length(inp)
    layout.(inp{ndx}) = inp{ndx+1};
end
layout.icon = layoutIcon(inp(3:2:length(inp)));


function icon = layoutIcon(arrList)

num = 128;
icon = ones(num,num);
for ndx = 1:length(arrList)
    pos = arrList{ndx};
    pos(3:4) = pos(1:2)+pos(3:4);
    pos(pos>1) = 1;
    pos = round(pos*(num-1))+1;
    [x y] = meshgrid(-1:2/(pos(3)-pos(1)):1,-1:2/(pos(4)-pos(2)):1);
    z = max(abs(x),abs(y)).^4;
    z = z/max(z(:));
    icon(pos(2):pos(4),pos(1):pos(3)) = z;
end
icon = imresize(icon,[16 16]);
icon = flipud((icon-min(icon(:)))/range(icon(:)));
icon = 1-cat(3,icon,icon,icon);