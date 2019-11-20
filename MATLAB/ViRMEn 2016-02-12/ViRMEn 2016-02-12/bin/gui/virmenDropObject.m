function virmenDropObject(varargin)

global virmenDragging

if isempty(virmenDragging) || ~isstruct(virmenDragging)
    virmenDragging = -1;
    return
end

global guifig;

set(guifig,'units','pixels');
pos = get(guifig,'currentpoint');
set(guifig,'units',virmenDragging.backupUnitsFigure);

df = pos - virmenDragging.startPt;

if norm(df)<2
    virmenDragging = [];
    return
end

set(virmenDragging.axes,'units','pixels');
axPos = get(virmenDragging.axes,'position');
set(virmenDragging.axes,'units',virmenDragging.backupUnitsAxes)

perc = df./axPos(3:4);
lm(1) = range(get(virmenDragging.axes,'xlim'));
lm(2) = range(get(virmenDragging.axes,'ylim'));
move = perc.*lm;

switch virmenDragging.type
    case 'shape'
        loc = virmenDragging.object.locations;
        loc = loc+repmat(move,size(loc,1),1);
        delete(virmenDragging.tempShape);
        virmenEventHandler('changeShapeLocations',loc);
    case 'object'
        loc = virmenDragging.object.locations;
        loc = loc+repmat(move,size(loc,1),1);
        delete(virmenDragging.tempShape);
        virmenEventHandler('changeObjectLocations',loc);
    case 'objectLocation'
        virmenDragging.x{virmenDragging.indx} = num2str(virmenDragging.startX(virmenDragging.indx) + move(1));
        virmenDragging.y{virmenDragging.indx} = num2str(virmenDragging.startY(virmenDragging.indx) + move(2));
        delete(virmenDragging.tempShape);
        virmenEventHandler('changeObjectLocations',[virmenDragging.x virmenDragging.y]);
    case 'shapeLocation'
        virmenDragging.x{virmenDragging.indx} = num2str(virmenDragging.startX(virmenDragging.indx) + move(1));
        virmenDragging.y{virmenDragging.indx} = num2str(virmenDragging.startY(virmenDragging.indx) + move(2));
        delete(virmenDragging.tempShape);
        virmenEventHandler('changeShapeLocations',[virmenDragging.x virmenDragging.y]);
    case 'shapeCopy'
        virmenDragging.startX(end+1,1) = virmenDragging.startX(virmenDragging.indx) + move(1);
        virmenDragging.startY(end+1,1) = virmenDragging.startY(virmenDragging.indx) + move(2);
        sz = length(virmenDragging.startX);
        if virmenDragging.indx == 1
            ord = [sz 1:sz-1];
        elseif virmenDragging.indx == length(virmenDragging.startX)
            ord = 1:sz;
        elseif norm([virmenDragging.startX(sz) virmenDragging.startY(sz)]-[virmenDragging.startX(virmenDragging.indx-1) virmenDragging.startY(virmenDragging.indx-1)]) < ...
                norm([virmenDragging.startX(sz) virmenDragging.startY(sz)]-[virmenDragging.startX(virmenDragging.indx+1) virmenDragging.startY(virmenDragging.indx+1)])
            ord = [1:virmenDragging.indx-1 sz virmenDragging.indx:sz-1];
        else
            ord = [1:virmenDragging.indx sz virmenDragging.indx+1:sz-1];
        end
        
        newX = virmenDragging.x;
        newY = virmenDragging.y;
        newX{end+1,1} = num2str(virmenDragging.startX(virmenDragging.indx) + move(1));
        newY{end+1,1} = num2str(virmenDragging.startY(virmenDragging.indx) + move(2));

        delete(virmenDragging.tempShape);
        virmenEventHandler('changeShapeLocations',[newX(ord) newY(ord)]);
    case 'objectCopy'
        virmenDragging.startX(end+1,1) = virmenDragging.startX(virmenDragging.indx) + move(1);
        virmenDragging.startY(end+1,1) = virmenDragging.startY(virmenDragging.indx) + move(2);
        sz = length(virmenDragging.startX);
        if virmenDragging.indx == 1
            ord = [sz 1:sz-1];
        elseif virmenDragging.indx == length(virmenDragging.startX)
            ord = 1:sz;
        elseif norm([virmenDragging.startX(sz) virmenDragging.startY(sz)]-[virmenDragging.startX(virmenDragging.indx-1) virmenDragging.startY(virmenDragging.indx-1)]) < ...
                norm([virmenDragging.startX(sz) virmenDragging.startY(sz)]-[virmenDragging.startX(virmenDragging.indx+1) virmenDragging.startY(virmenDragging.indx+1)])
            ord = [1:virmenDragging.indx-1 sz virmenDragging.indx:sz-1];
        else
            ord = [1:virmenDragging.indx sz virmenDragging.indx+1:sz-1];
        end
        
        newX = virmenDragging.x;
        newY = virmenDragging.y;
        newX{end+1,1} = num2str(virmenDragging.startX(virmenDragging.indx) + move(1));
        newY{end+1,1} = num2str(virmenDragging.startY(virmenDragging.indx) + move(2));
        
        delete(virmenDragging.tempShape);
        virmenEventHandler('changeObjectLocations',[newX(ord) newY(ord)]);
end

virmenDragging = [];