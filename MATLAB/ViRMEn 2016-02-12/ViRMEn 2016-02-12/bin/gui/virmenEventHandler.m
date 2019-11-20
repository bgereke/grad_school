function virmenEventHandler(type,inp)

global virmenDragging
virmenDragging = [];

global guifig; 
currax = get(guifig,'currentaxes');

switch type
    case 'selectShape'
        %
    case 'selectObject'
        %
    case 'zoomOnWorld'
        %
    case 'zoomOnTexture'
        %
    otherwise
        set(guifig,'Pointer','watch');
        drawnow
end

handles = guidata(guifig);
if isempty(handles)
    return
end

wNum = handles.state.selectedWorld;
oNum = handles.state.selectedObject;
sNum = handles.state.selectedShape;

justSaved = false;

handles.exper.updateCodeText;

if handles.historyBool.(type)
    undo.param = {};
    redo.param = {};
    oldState = handles.state;
end

switch type
    case 'changeView'
        ax = findobj(handles.figs.worldDrawing,'type','axes');
        if ~isempty(ax)
            switch inp
                case 'isometric'
                    set(ax,'view',[-45 45]);
                case 'top'
                    set(ax,'view',[0 90]);
                case 'front'
                    set(ax,'view',[0 0]);
                case 'side'
                    set(ax,'view',[90 0]);
            end
        end
    case 'rotateView';
        ax = findobj(handles.figs.worldDrawing,'type','axes');
        if ~isempty(ax)
            v = get(ax,'view');
            v = round(v/(45/3))*(45/3);
            switch inp
                case 'right'
                    v(1) = v(1)-45/3;
                case 'left'
                    v(1) = v(1)+45/3;
                case 'down'
                    v(2) = v(2)+45/3;
                    if v(2)>90
                        v(2) = 90;
                    end
                case 'up'
                    v(2) = v(2)-45/3;
                    if v(2)<-90
                        v(2) = -90;
                    end
            end
            set(ax,'view',v);
        end
    case 'editCode'
        if strcmp(func2str(handles.exper.(inp)),'undefined')
            errordlg('Cannot edit undefined code.','Error');
        else
            if exist([func2str(handles.exper.(inp)) '.m'],'file')
                fname = [func2str(handles.exper.(inp)) '.m'];
            elseif exist([func2str(handles.exper.(inp)) '.c'],'file')
                fname = [func2str(handles.exper.(inp)) '.c'];
            else
                fname = [];
            end
            if ~isempty(fname)
                edit(fname);
            else
                errordlg(['Cannot edit ' func2str(handles.exper.(inp)) '.'],'Error');
            end
        end
        set(guifig,'Pointer','arrow'); return
    case 'runMacro'
        try
            if ~isfield(handles.macros,inp)
                handles.macros.(inp) = eval(inp);
            end
            mcr = handles.macros.(inp);
        
            figs = fieldnames(handles.figNames);
            for ndx = 1:length(figs)
                handles.bools(handles.evtNames.('runMacro'),handles.figNames.(figs{ndx})) = mcr.(['update' upper(figs{ndx}(1)) figs{ndx}(2:end)]);
            end
            
            if mcr.trackHistory
                undo.param = {'runMacro',handles.exper,false};
                exper = copyVirmenObject(handles.exper);
                mcr.exper = exper;
                mcr.guiState = handles.state;
                mcr.run;
                handles.state = mcr.guiState;
                handles.exper = exper;
                redo.param = {'runMacro',handles.exper,false};
            else
                undo.param = {'runMacro',handles.exper,true};
                mcr.exper = handles.exper;
                mcr.guiState = handles.state;
                mcr.run;
                handles.state = mcr.guiState;
                redo.param = {'runMacro',handles.exper,true};
            end
            
        catch err
            errordlg(['Macro function ' err.stack(1).name ' generated an error on line ' num2str(err.stack(1).line) ': ' err.message],'Error');
            set(guifig,'Pointer','arrow'); return
        end

    case 'changeMacroSettings'
        try
            if ~isfield(handles.macros,inp)
                handles.macros.(inp) = eval(inp);
            end
            mcr = handles.macros.(inp);
            
            props = setdiff(properties(mcr),properties(virmenMacro));
            propVal1 = struct;
            for ndx = 1:length(props)
                propVal1.(props{ndx}) = mcr.(props{ndx});
            end
            
            mcr.exper = handles.exper;
            mcr.guiState = handles.state;
            mcr.settings;
            
            propVal2 = struct;
            for ndx = 1:length(props)
                propVal2.(props{ndx}) = mcr.(props{ndx});
            end
            if ~isequal(propVal1,propVal2)
                undo.param = {'macroSettings',inp,propVal1};
                redo.param = {'macroSettings',inp,propVal2};
            else
                set(guifig,'Pointer','arrow'); return
            end
        catch err
            errordlg(['Macro function ' err.stack(1).name ' generated an error on line ' num2str(err.stack(1).line) ': ' err.message],'Error');
            set(guifig,'Pointer','arrow'); return
        end
        
    case 'startProgram'
        [handles.figs, handles.buttons, handles.menus] = createFigures;
        handles.separated = setdiff(findall(guifig,'separator','on'),findall(guifig,'type','uimenu'));
        justSaved = true;
        
        undo.param = {'startProgram'};
        redo.param = {'startProgram'};
        
        set(guifig,'resizefcn','virmenEventHandler(''resizeFigure'',''n/a'')');
    case 'selectShape'
        if ischar(inp)
            inp = {sNum, 'none', 1};
            str = 'color';
        else
            str = get(guifig,'selectiontype');
        end
        switch str
            case 'normal'
                handles.state.selectedShape = inp{1};
                if handles.state.selectedShape > 1
                    if isempty(virmenDragging)
                        virmenDragging.type = 'shape';
                        virmenDragging.object = handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{inp{1}};
                        virmenDragging.backupUnitsFigure = get(guifig,'units');
                        virmenDragging.backupUnitsAxes = get(currax,'units');
                        set(guifig,'units','pixels');
                        virmenDragging.startPt = get(guifig,'currentpoint');
                        set(guifig,'units',virmenDragging.backupUnitsFigure);
                        [x, y] = handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{inp{1}}.coords2D;
                        virmenDragging.axes = currax;
                        virmenDragging.x = x;
                        virmenDragging.y = y;
                        virmenDragging.tempShape = [];
                        virmenDragging.marker = inp{2};
                        virmenDragging.markerSize = inp{3};
                    else
                        virmenDragging = [];
                    end
                end
            case 'extend'
                if handles.state.selectedShape > 1
                    virmenDragging.type = 'shapeLocation';
                    virmenDragging.object = copyVirmenObject(handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{inp{1}});
                    virmenDragging.backupUnitsFigure = get(guifig,'units');
                    virmenDragging.backupUnitsAxes = get(currax,'units');
                    pos = get(currax,'currentpoint');
                    pos = pos(1,1:2);
                    set(guifig,'units','pixels');
                    virmenDragging.startPt = get(guifig,'currentpoint');
                    set(guifig,'units',virmenDragging.backupUnitsFigure);
                    virmenDragging.axes = currax;
                    loc = handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{inp{1}}.locations;
                    virmenDragging.startX = loc(:,1);
                    virmenDragging.startY = loc(:,2);
                    virmenDragging.x = handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{inp{1}}.symbolic.x;
                    virmenDragging.y = handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{inp{1}}.symbolic.y;
                    if ~iscell(virmenDragging.x)
                        virmenDragging.x = cell(size(loc,1),1);
                        for ndx = 1:size(loc,1)
                            virmenDragging.x{ndx} = num2str(loc(ndx,1));
                        end
                    end
                    if ~iscell(virmenDragging.y)
                        virmenDragging.y = cell(size(loc,1),1);
                        for ndx = 1:size(loc,1)
                            virmenDragging.y{ndx} = num2str(loc(ndx,2));
                        end
                    end
                    virmenDragging.tempShape = [];
                    virmenDragging.marker = inp{2};
                    virmenDragging.markerSize = inp{3};
                    dst = sum(bsxfun(@minus,loc,pos).^2,2);
                    [~, virmenDragging.indx] = min(dst);
                end
            case 'alt'
                if handles.state.selectedShape > 1
                    virmenDragging.type = 'shapeCopy';
                    virmenDragging.object = copyVirmenObject(handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{inp{1}});
                    virmenDragging.backupUnitsFigure = get(guifig,'units');
                    virmenDragging.backupUnitsAxes = get(currax,'units');
                    pos = get(currax,'currentpoint');
                    pos = pos(1,1:2);
                    set(guifig,'units','pixels');
                    virmenDragging.startPt = get(guifig,'currentpoint');
                    set(guifig,'units',virmenDragging.backupUnitsFigure);
                    virmenDragging.axes = currax;
                    loc = handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{inp{1}}.locations;
                    virmenDragging.startX = loc(:,1);
                    virmenDragging.startY = loc(:,2);
                    virmenDragging.x = handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{inp{1}}.symbolic.x;
                    virmenDragging.y = handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{inp{1}}.symbolic.y;
                    if ~iscell(virmenDragging.x)
                        virmenDragging.x = cell(size(loc,1),1);
                        for ndx = 1:size(loc,1)
                            virmenDragging.x{ndx} = num2str(loc(ndx,1));
                        end
                    end
                    if ~iscell(virmenDragging.y)
                        virmenDragging.y = cell(size(loc,1),1);
                        for ndx = 1:size(loc,1)
                            virmenDragging.y{ndx} = num2str(loc(ndx,2));
                        end
                    end
                    virmenDragging.tempShape = [];
                    virmenDragging.marker = inp{2};
                    virmenDragging.markerSize = inp{3};
                    dst = sum(bsxfun(@minus,loc,pos).^2,2);
                    [~, virmenDragging.indx] = min(dst);
                end
            case {'open','color'}
                undo.param = {'GUIOperation',handles.state};
                handles.state.selectedShape = inp{1};
                redo.param = {'GUIOperation',handles.state};
                if oNum == 0
                    errordlg('No color selected.','Error');
                    set(guifig,'Pointer','arrow'); return
                end
                if isa(handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{inp{1}},'shapeColor')
                    col = handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{inp{1}}.RGBA;
                    newcol = uisetcolor(col(1:3));
                    str = {num2str(newcol(1)); num2str(newcol(2)); num2str(newcol(3)); num2str(col(4))};
                    guidata(guifig, handles)
                    virmenEventHandler('changeShapeProperties',{{'R','G','B','Alpha'},str});
                    set(guifig,'Pointer','arrow'); return
                elseif strcmp(str,'color')
                    errordlg('No color selected.','Error');
                    set(guifig,'Pointer','arrow'); return
                end
        end
        
    case 'addShape'
        if oNum == 0
            errordlg('No object with texture selected.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        if length(handles.exper.worlds{wNum}.objects{oNum}.texture.shapes)==2 ...
                & isa(handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{2},'shapeImage') %#ok<AND2>
            errordlg('Cannot add shapes to an image texture.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        subplot(findobj(handles.figs.textureSketch,'type','axes'))
        obj = eval(inp);
        try
            set(guifig,'handlevisibility','on');
            set(get(gca,'title'),'string',obj.helpString,'fontsize',10);
            obj = getPoints(obj);
            set(get(gca,'title'),'string','');
            set(guifig,'handlevisibility','off');
        catch %#ok<CTCH>
            errordlg('Error entering shape locations.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        if ~isempty(obj.locations)
            handles.exper.worlds{wNum}.objects{oNum}.texture = addShape(handles.exper.worlds{wNum}.objects{oNum}.texture,obj,'new');
            handles.state.selectedShape = length(handles.exper.worlds{wNum}.objects{oNum}.texture.shapes);
            undo.param = {'deleteShape',wNum,oNum,handles.state.selectedShape};
            redo.param = {'addShape',wNum,oNum,handles.state.selectedShape,copyVirmenObject(obj)};
        else
            set(guifig,'Pointer','arrow'); return
        end
        
    case 'deleteShape'
        if oNum == 0
            errordlg('No shape selected.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        if sNum == 1
            errordlg('Cannot delete boundary.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        undo.param = {'addShape',wNum,oNum,sNum, handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}};
        redo.param = {'deleteShape',wNum,oNum,sNum};
        handles.exper.worlds{wNum}.objects{oNum}.texture.shapes(sNum) = [];
        if handles.state.selectedShape > length(handles.exper.worlds{wNum}.objects{oNum}.texture.shapes)
            handles.state.selectedShape = handles.state.selectedShape - 1;
        end
    case 'computeTexture'
        if oNum == 0
            errordlg('No object with texture selected.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        oth = zeros(0,2);
        for w = 1:length(handles.exper.worlds)
            for o = 1:length(handles.exper.worlds{w}.objects)
                if ~all([w o]==[wNum oNum]) && isequalwithequalnans(handles.exper.worlds{w}.objects{o}.texture.triangles, ...
                        handles.exper.worlds{wNum}.objects{oNum}.texture.triangles) %#ok<DISEQN>
                    oth(end+1,:) = [w o]; %#ok<AGROW>
                end
            end
        end
        undo.param = {{'changeProperty',[wNum oNum],'texture',copyVirmenObject(handles.exper.worlds{wNum}.objects{oNum}.texture)}};
        handles.exper.worlds{wNum}.objects{oNum}.texture.compute;
        redo.param = {{'changeProperty',[wNum oNum],'texture',copyVirmenObject(handles.exper.worlds{wNum}.objects{oNum}.texture)}};
        if ~isempty(oth)
            str = cell(1,size(oth,1));
            for ndx = 1:size(oth,1)
                str{ndx} = handles.exper.worlds{oth(ndx,1)}.objects{oth(ndx,2)}.fullName;
            end
            [indx,ok] = listdlg('ListString',str,'InitialValue',1:length(str),'ListSize',[250 150], ...
                'Name','Select objects','PromptString','Change other objects with identical texture:');
            if ok > 0
                for ndx = 1:length(indx)
                    undo.param{end+1,1} = {'changeProperty',oth(indx(ndx),:),'texture',copyVirmenObject(handles.exper.worlds{oth(indx(ndx),1)}.objects{oth(indx(ndx),2)}.texture)};
                    setTexture(handles.exper.worlds{oth(indx(ndx),1)}.objects{oth(indx(ndx),2)},handles.exper.worlds{wNum}.objects{oNum}.texture,'copy');
                    redo.param{end+1,1} = {'changeProperty',oth(indx(ndx),:),'texture',copyVirmenObject(handles.exper.worlds{oth(indx(ndx),1)}.objects{oth(indx(ndx),2)}.texture)};
                end
            end
        end
        if size(undo.param,1)==1
            undo.param = undo.param{1};
            redo.param = redo.param{1};
        end
    case 'refineTexture'
        if oNum == 0
            errordlg('No object with texture selected.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        if length(handles.exper.worlds{wNum}.objects{oNum}.texture.shapes)==2 ...
                & isa(handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{2},'shapeImage') %#ok<AND2>
            errordlg('To refine this texture, change image size.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        answer = refiningGui(handles.exper.worlds{wNum}.objects{oNum}.texture,handles.state.triangulationColor);
        if isempty(answer)
            set(guifig,'Pointer','arrow'); return
        end
        virmenEventHandler('changeTextureRefining',answer);
        set(guifig,'Pointer','arrow'); return
    case 'changeTextureRefining'
        oth = zeros(0,2);
        for w = 1:length(handles.exper.worlds)
            for o = 1:length(handles.exper.worlds{w}.objects)
                if ~all([w o]==[wNum oNum]) && isequalwithequalnans(handles.exper.worlds{w}.objects{o}.texture.triangles, ...
                        handles.exper.worlds{wNum}.objects{oNum}.texture.triangles) %#ok<DISEQN>
                    oth(end+1,:) = [w o]; %#ok<AGROW>
                end
            end
        end
        undo.param = {{'changeProperty',[wNum oNum],'texture',copyVirmenObject(handles.exper.worlds{wNum}.objects{oNum}.texture)}};
        handles.exper.worlds{wNum}.objects{oNum}.texture.refining = inp(1:2)';
        handles.exper.worlds{wNum}.objects{oNum}.texture.grid = inp(3:4)';
        handles.exper.worlds{wNum}.objects{oNum}.texture.tilable = inp(5:6)';
        handles.exper.worlds{wNum}.objects{oNum}.texture.compute;
        redo.param = {{'changeProperty',[wNum oNum],'texture',copyVirmenObject(handles.exper.worlds{wNum}.objects{oNum}.texture)}};
        if ~isempty(oth)
            str = cell(1,size(oth,1));
            for ndx = 1:size(oth,1)
                str{ndx} = handles.exper.worlds{oth(ndx,1)}.objects{oth(ndx,2)}.fullName;
            end
            [indx,ok] = listdlg('ListString',str,'InitialValue',1:length(str),'ListSize',[250 150], ...
                'Name','Select objects','PromptString','Change other objects with identical texture:');
            if ok > 0
                for ndx = 1:length(indx)
                    undo.param{end+1,1} = {'changeProperty',oth(indx(ndx),:),'texture',copyVirmenObject(handles.exper.worlds{oth(indx(ndx),1)}.objects{oth(indx(ndx),2)}.texture)};
                    setTexture(handles.exper.worlds{oth(indx(ndx),1)}.objects{oth(indx(ndx),2)},handles.exper.worlds{wNum}.objects{oNum}.texture,'copy');
                    redo.param{end+1,1} = {'changeProperty',oth(indx(ndx),:),'texture',copyVirmenObject(handles.exper.worlds{oth(indx(ndx),1)}.objects{oth(indx(ndx),2)}.texture)};
                end
            end
        end
        if size(undo.param,1)==1
            undo.param = undo.param{1};
            redo.param = redo.param{1};
        end
    case 'changeTriangulationVisibility'
        if strcmp(inp{1},'on')
            handles.state.showTriangulation = 1;
        elseif strcmp(inp{1},'off')
            handles.state.showTriangulation = 0;
        else
            handles.state.showTriangulation = 1-handles.state.showTriangulation;
        end
    case 'changeTriangulationColor'
        handles.state.triangulationColor = uisetcolor(handles.state.triangulationColor);
    case 'zoomOnTexture'
        switch get(guifig,'selectiontype')
            case 'normal'
                set(guifig,'Pointer','custom','PointerShapeCdata',zoomPointer);
                rect = virmenZoom(guifig,currax);
                set(guifig,'Pointer','watch');
                
                xl = [rect(1) rect(1)+rect(3)];
                yl = [rect(2) rect(2)+rect(4)];

                if rect(3)==0 || rect(4)==0
                    set(guifig,'Pointer','arrow'); return
                end
            case 'extend'
                buff = 1.05;
                xl = [0 handles.exper.worlds{wNum}.objects{oNum}.texture.width];
                yl = [0 handles.exper.worlds{wNum}.objects{oNum}.texture.height];
                xl = [mean(xl)-buff/2*range(xl) mean(xl)+buff/2*range(xl)];
                yl = [mean(yl)-buff/2*range(yl) mean(yl)+buff/2*range(yl)];
            case 'alt'
                pos = get(currax,'currentpoint');
                pos = pos(1,1:2);
                xl = get(currax,'xlim');
                yl = get(currax,'ylim');
                buff = 2;
                xl = [pos(1)-buff/2*range(xl) pos(1)+buff/2*range(xl)];
                yl = [pos(2)-buff/2*range(yl) pos(2)+buff/2*range(yl)];
            otherwise
                set(guifig,'Pointer','arrow'); return
        end
        
        handles.state.textureXLim = xl;
        handles.state.textureYLim = yl;
    case 'zoomOnWorld'
        switch get(guifig,'selectiontype')
            case 'normal'
                set(guifig,'Pointer','custom','PointerShapeCdata',zoomPointer);
                rect = virmenZoom(guifig,currax);
                set(guifig,'Pointer','watch');
                
                xl = [rect(1) rect(1)+rect(3)];
                yl = [rect(2) rect(2)+rect(4)];

                if rect(3)==0 || rect(4)==0
                    set(guifig,'Pointer','arrow'); return
                end
            case 'extend'
                if isempty(handles.exper.worlds{wNum}.objects)
                    xl = handles.defaultProperties.worldXLim;
                    yl = handles.defaultProperties.worldYLim;
                else
                    axis(currax,'tight');
                    xl = get(currax,'xlim');
                    yl = get(currax,'ylim');
                    buff = 1.05;
                    xl = [mean(xl)-buff/2*range(xl) mean(xl)+buff/2*range(xl)];
                    yl = [mean(yl)-buff/2*range(yl) mean(yl)+buff/2*range(yl)];
                end
            case 'alt'
                pos = get(currax,'currentpoint');
                pos = pos(1,1:2);
                xl = get(currax,'xlim');
                yl = get(currax,'ylim');
                buff = 2;
                xl = [pos(1)-buff/2*range(xl) pos(1)+buff/2*range(xl)];
                yl = [pos(2)-buff/2*range(yl) pos(2)+buff/2*range(yl)];
            otherwise
                set(guifig,'Pointer','arrow'); return
        end
        
        handles.state.worldXLim = xl;
        handles.state.worldYLim = yl;
        
    case 'selectObject'
        handles.state.selectedObject = inp{1};
        handles.state.selectedShape = 1;
        if strcmp(get(guifig,'selectiontype'),'normal') && handles.state.selectedObject > 0 && get(currax,'parent')==handles.worldSketch
            if isempty(virmenDragging)
                virmenDragging.type = 'object';
                virmenDragging.object = handles.exper.worlds{wNum}.objects{handles.state.selectedObject};
                virmenDragging.backupUnitsFigure = get(guifig,'units');
                virmenDragging.backupUnitsAxes = get(currax,'units');
                set(guifig,'units','pixels');
                virmenDragging.startPt = get(guifig,'currentpoint');
                set(guifig,'units',virmenDragging.backupUnitsFigure);
                [x, y, ~] = handles.exper.worlds{wNum}.objects{handles.state.selectedObject}.coords2D;
                virmenDragging.axes = currax;
                virmenDragging.x = x;
                virmenDragging.y = y;
                virmenDragging.tempShape = [];
                virmenDragging.marker = inp{2};
                virmenDragging.markerSize = inp{3};
            else
                virmenDragging = [];
            end
        elseif strcmp(get(guifig,'selectiontype'),'extend') && handles.state.selectedObject > 0 && get(currax,'parent')==handles.worldSketch
            virmenDragging.type = 'objectLocation';
            virmenDragging.object = copyVirmenObject(handles.exper.worlds{wNum}.objects{handles.state.selectedObject});
            virmenDragging.backupUnitsFigure = get(guifig,'units');
            virmenDragging.backupUnitsAxes = get(currax,'units');
            pos = get(currax,'currentpoint');
            pos = pos(1,1:2);
            set(guifig,'units','pixels');
            virmenDragging.startPt = get(guifig,'currentpoint');
            set(guifig,'units',virmenDragging.backupUnitsFigure);
            virmenDragging.axes = currax;
            loc = handles.exper.worlds{wNum}.objects{handles.state.selectedObject}.locations;
            virmenDragging.startX = loc(:,1);
            virmenDragging.startY = loc(:,2);
            virmenDragging.x = handles.exper.worlds{wNum}.objects{handles.state.selectedObject}.symbolic.x;
            virmenDragging.y = handles.exper.worlds{wNum}.objects{handles.state.selectedObject}.symbolic.y;
            if ~iscell(virmenDragging.x)
                virmenDragging.x = cell(size(loc,1),1);
                for ndx = 1:size(loc,1)
                    virmenDragging.x{ndx} = num2str(loc(ndx,1));
                end
            end
            if ~iscell(virmenDragging.y)
                virmenDragging.y = cell(size(loc,1),1);
                for ndx = 1:size(loc,1)
                    virmenDragging.y{ndx} = num2str(loc(ndx,2));
                end
            end
            virmenDragging.tempShape = [];
            virmenDragging.marker = inp{2};
            virmenDragging.markerSize = inp{3};
            dst = sum(bsxfun(@minus,loc,pos).^2,2);
            [~, virmenDragging.indx] = min(dst);
        elseif strcmp(get(guifig,'selectiontype'),'alt') && handles.state.selectedObject > 0 && get(currax,'parent')==handles.worldSketch
            virmenDragging.type = 'objectCopy';
            virmenDragging.object = copyVirmenObject(handles.exper.worlds{wNum}.objects{handles.state.selectedObject});
            virmenDragging.backupUnitsFigure = get(guifig,'units');
            virmenDragging.backupUnitsAxes = get(currax,'units');
            pos = get(currax,'currentpoint');
            pos = pos(1,1:2);
            set(guifig,'units','pixels');
            virmenDragging.startPt = get(guifig,'currentpoint');
            set(guifig,'units',virmenDragging.backupUnitsFigure);
            virmenDragging.axes = currax;
            loc = handles.exper.worlds{wNum}.objects{handles.state.selectedObject}.locations;
            virmenDragging.startX = loc(:,1);
            virmenDragging.startY = loc(:,2);
            virmenDragging.x = handles.exper.worlds{wNum}.objects{handles.state.selectedObject}.symbolic.x;
            virmenDragging.y = handles.exper.worlds{wNum}.objects{handles.state.selectedObject}.symbolic.y;
            if ~iscell(virmenDragging.x)
                virmenDragging.x = cell(size(loc,1),1);
                for ndx = 1:size(loc,1)
                    virmenDragging.x{ndx} = num2str(loc(ndx,1));
                end
            end
            if ~iscell(virmenDragging.y)
                virmenDragging.y = cell(size(loc,1),1);
                for ndx = 1:size(loc,1)
                    virmenDragging.y{ndx} = num2str(loc(ndx,2));
                end
            end
            virmenDragging.tempShape = [];
            virmenDragging.marker = inp{2};
            virmenDragging.markerSize = inp{3};
            dst = sum(bsxfun(@minus,loc,pos).^2,2);
            [~, virmenDragging.indx] = min(dst);
        elseif strcmp(get(guifig,'selectiontype'),'open') && handles.state.selectedObject > 0
            virmenEventHandler('builtinLayout','texture');
            set(guifig,'Pointer','arrow'); return
        end
        
    case 'addObject'
        subplot(findobj(handles.figs.worldSketch,'type','axes'))
        obj = eval(inp);
        try
            set(guifig,'handlevisibility','on');
            set(get(gca,'title'),'string',obj.helpString,'fontsize',10);
            obj = getPoints(obj);
            set(get(gca,'title'),'string','');
            set(guifig,'handlevisibility','off');
        catch %#ok<CTCH>
            errordlg('Error entering object locations.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        obj.tiling = handles.defaultProperties.tiling;
        obj.edgeRadius = handles.defaultProperties.edgeRadius;
        obj.texture.tilable = handles.defaultProperties.textureTilable;
        obj.texture.refining = handles.defaultProperties.triangulationRefining;
        obj.texture.grid = handles.defaultProperties.triangulationGrid;
        obj.texture.compute;
        
        if ~isempty(obj.locations)
            handles.exper.worlds{wNum} = addObject(handles.exper.worlds{wNum},obj,'new');
            handles.state.selectedObject = length(handles.exper.worlds{wNum}.objects);
            handles.state.selectedShape = 1;
            undo.param = {'deleteObject',wNum,handles.state.selectedObject};
            redo.param = {'addObject',wNum,handles.state.selectedObject,copyVirmenObject(obj)};
        else
            set(guifig,'Pointer','arrow'); return
        end
    case 'deleteObject'
        if oNum == 0
            errordlg('Cannot delete initial conditions indicator.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        undo.param = {'addObject',wNum,oNum,handles.exper.worlds{wNum}.objects{oNum}};
        redo.param = {'deleteObject',wNum,oNum};
        handles.exper.worlds{wNum}.objects(oNum) = [];
        if handles.state.selectedObject > length(handles.exper.worlds{wNum}.objects)
            handles.state.selectedObject = handles.state.selectedObject - 1;
        end
    case 'changeWorldBackground'
        undo.param = {'changeProperty',wNum,'backgroundColor',handles.exper.worlds{wNum}.getValue.backgroundColor};
        c = uisetcolor(handles.exper.worlds{wNum}.backgroundColor);
        handles.exper.worlds{wNum}.backgroundColor = c;
        redo.param = {'changeProperty',wNum,'backgroundColor',handles.exper.worlds{wNum}.getValue.backgroundColor};
    case 'changeWorldTransparency'
        undo.param = {'changeProperty',wNum,'transparency',handles.exper.worlds{wNum}.getValue.transparency};
        switch inp{1}
            case 'on'
                handles.exper.worlds{wNum}.transparency = 1;
            case 'off'
                handles.exper.worlds{wNum}.transparency = 0;
            case 'switch'
                handles.exper.worlds{wNum}.transparency = 1-handles.exper.worlds{wNum}.transparency;
                
        end
        redo.param = {'changeProperty',wNum,'transparency',handles.exper.worlds{wNum}.getValue.transparency};
    case 'switchWireframe'
        switch inp{1}
            case 'on'
                handles.state.showWireframe = 1;
            case 'off'
                handles.state.showWireframe = 0;
            case 'switch'
                handles.state.showWireframe = 1-handles.state.showWireframe;
                
        end
    case 'changeTiling'
        undo.param = {'changeProperty',[wNum oNum],'tiling',handles.exper.worlds{wNum}.objects{oNum}.getValue.tiling};
        til = handles.exper.worlds{wNum}.objects{oNum}.getValue.tiling;
        til{inp{1}} = inp{2};
        handles.exper.worlds{wNum}.objects{oNum}.tiling = til;
        redo.param = {'changeProperty',[wNum oNum],'tiling',handles.exper.worlds{wNum}.objects{oNum}.getValue.tiling};
    case 'renameObject'
        if oNum == 0
            errordlg('Cannot rename initial conditions indicator.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        undo.param = {{'changeProperty',[wNum oNum],'name',handles.exper.worlds{wNum}.objects{oNum}.name};
            {'changeProperty',[wNum oNum NaN],'name',handles.exper.worlds{wNum}.objects{oNum}.texture.name}};
        answer = inputdlg({['New name for object ' handles.exper.worlds{wNum}.objects{oNum}.fullName],...
            'New name for the object''s texture'}, ...
            'New name',1,{handles.exper.worlds{wNum}.objects{oNum}.name, handles.exper.worlds{wNum}.objects{oNum}.texture.name});
        if isempty(answer)
            set(guifig,'Pointer','arrow'); return
        end
        if ~isvarname(answer{1})
            errordlg(['''' answer{1} ''' is an invalid variable name.'],'Error','modal');
            set(guifig,'Pointer','arrow'); return
        end
        if ~isvarname(answer{2})
            errordlg(['''' answer{2} ''' is an invalid variable name.'],'Error','modal');
            set(guifig,'Pointer','arrow'); return
        end
        handles.exper.worlds{wNum}.objects{oNum}.name = answer{1};
        handles.exper.worlds{wNum}.objects{oNum}.texture.name = answer{2};
        redo.param = {{'changeProperty',[wNum oNum],'name',handles.exper.worlds{wNum}.objects{oNum}.name};
            {'changeProperty',[wNum oNum NaN],'name',handles.exper.worlds{wNum}.objects{oNum}.texture.name}};
    case 'sortObjects'
        if length(handles.exper.worlds{wNum}.objects)<2
            errordlg('Need at least two objects to reorder.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        fld = cell(1,length(handles.exper.worlds{wNum}.objects));
        for ndx = 1:length(fld)
            fld{ndx} = handles.exper.worlds{wNum}.objects{ndx}.fullName;
        end
        ord = reorderVariables(fld,'Objects');
        if ~all(ord==(1:length(ord)))
            [~,orig] = sort(ord);
            undo.param = {'sortObjects',wNum,orig};
            redo.param = {'sortObjects',wNum,ord};
            handles.exper.worlds{wNum}.objects = handles.exper.worlds{wNum}.objects(ord);
            if oNum > 0
                handles.state.selectedObject = find(ord==oNum);
            end
        else
            set(guifig,'Pointer','arrow'); return
        end
    case 'changeObjectProperties'
        if oNum == 0
            switch inp{1}
                case {'X','Y','Z','rotation'}
                    str = handles.exper.worlds{wNum}.getValue.startLocation;
                    switch inp{1}
                        case 'X'
                            str{1} = inp{2};
                        case 'Y'
                            str{2} = inp{2};
                        case 'Z'
                            str{3} = inp{2};
                        case 'rotation'
                            str{4} = inp{2};
                    end
                    undo.param = {'changeProperty',wNum,'startLocation',handles.exper.worlds{wNum}.getValue.startLocation};
                    handles.exper.worlds{wNum}.startLocation = str;
                    redo.param = {'changeProperty',wNum,'startLocation',handles.exper.worlds{wNum}.getValue.startLocation};
                case {'backgroundR','backgroundG','backgroundB'}
                    str = handles.exper.worlds{wNum}.getValue.backgroundColor;
                    switch inp{1}
                        case 'backgroundR'
                            str{1} = inp{2};
                        case 'backgroundG'
                            str{2} = inp{2};
                        case 'backgroundB'
                            str{3} = inp{2};
                    end
                    undo.param = {'changeProperty',wNum,'backgroundColor',handles.exper.worlds{wNum}.getValue.backgroundColor};
                    handles.exper.worlds{wNum}.backgroundColor = str;
                    redo.param = {'changeProperty',wNum,'backgroundColor',handles.exper.worlds{wNum}.getValue.backgroundColor};
                case 'enableTransparency'
                    undo.param = {'changeProperty',wNum,'transparency',handles.exper.worlds{wNum}.getValue.transparency};
                    handles.exper.worlds{wNum}.transparency = inp{2};
                    redo.param = {'changeProperty',wNum,'transparency',handles.exper.worlds{wNum}.getValue.transparency};
            end
        else
            undo.param = {'changeProperty',[wNum oNum],inp{1},handles.exper.worlds{wNum}.objects{oNum}.getValue.(inp{1})};
            handles.exper.worlds{wNum}.objects{oNum}.(inp{1}) = inp{2};
            redo.param = {'changeProperty',[wNum oNum],inp{1},handles.exper.worlds{wNum}.objects{oNum}.getValue.(inp{1})};
        end
    case 'changeObjectLocations'
        if ischar(handles.exper.worlds{wNum}.objects{oNum}.symbolic.x)
            undo.param = {{'changeProperty',[wNum oNum],'x',handles.exper.worlds{wNum}.objects{oNum}.symbolic.x};
                {'changeProperty',[wNum oNum],'y',handles.exper.worlds{wNum}.objects{oNum}.symbolic.y}};
        else
            undo.param = {{'changeProperty',[wNum oNum],'x',handles.exper.worlds{wNum}.objects{oNum}.getValue.x};
                {'changeProperty',[wNum oNum],'y',handles.exper.worlds{wNum}.objects{oNum}.getValue.y}};
        end
        handles.exper.worlds{wNum}.objects{oNum}.x = inp(:,1);
        handles.exper.worlds{wNum}.objects{oNum}.y = inp(:,2);
        redo.param = {{'changeProperty',[wNum oNum],'x',handles.exper.worlds{wNum}.objects{oNum}.getValue.x};
            {'changeProperty',[wNum oNum],'y',handles.exper.worlds{wNum}.objects{oNum}.getValue.y}};
    case 'addObjectLocation'
        if oNum == 0
            errordlg('Use properties to motify initial conditions.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        val = get(handles.table_objectLocations,'data');
        val = val(:,2:end);
        if ischar(handles.exper.worlds{wNum}.objects{oNum}.symbolic.x)
            undo.param = {{'changeProperty',[wNum oNum],'x',handles.exper.worlds{wNum}.objects{oNum}.symbolic.x};
                {'changeProperty',[wNum oNum],'y',handles.exper.worlds{wNum}.objects{oNum}.symbolic.y}};
        else
            undo.param = {{'changeProperty',[wNum oNum],'x',handles.exper.worlds{wNum}.objects{oNum}.getValue.x};
                {'changeProperty',[wNum oNum],'y',handles.exper.worlds{wNum}.objects{oNum}.getValue.y}};
        end
        val = [val; {'0','0'}];
        handles.exper.worlds{wNum}.objects{oNum}.x = val(:,1);
        handles.exper.worlds{wNum}.objects{oNum}.y = val(:,2);
        redo.param = {{'changeProperty',[wNum oNum],'x',handles.exper.worlds{wNum}.objects{oNum}.getValue.x};
            {'changeProperty',[wNum oNum],'y',handles.exper.worlds{wNum}.objects{oNum}.getValue.y}};
    case 'deleteObjectLocations'
        if oNum == 0
            errordlg('Use properties to motify initial conditions.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        indx = get(handles.table_objectLocations,'userdata');
        if isempty(indx)
            errordlg('No locations selected.','Error')
            set(guifig,'Pointer','arrow'); return
        end
        val = get(handles.table_objectLocations,'data');
        val = val(:,2:end);
        if size(val,1) == 1
            errordlg('At least one location is required.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        if ischar(handles.exper.worlds{wNum}.objects{oNum}.symbolic.x)
            undo.param = {{'changeProperty',[wNum oNum],'x',handles.exper.worlds{wNum}.objects{oNum}.symbolic.x};
                {'changeProperty',[wNum oNum],'y',handles.exper.worlds{wNum}.objects{oNum}.symbolic.y}};
        else
            undo.param = {{'changeProperty',[wNum oNum],'x',handles.exper.worlds{wNum}.objects{oNum}.getValue.x};
                {'changeProperty',[wNum oNum],'y',handles.exper.worlds{wNum}.objects{oNum}.getValue.y}};
        end
        val(indx,:) = [];
        handles.exper.worlds{wNum}.objects{oNum}.x = val(:,1);
        handles.exper.worlds{wNum}.objects{oNum}.y = val(:,2);
        redo.param = {{'changeProperty',[wNum oNum],'x',handles.exper.worlds{wNum}.objects{oNum}.getValue.x};
            {'changeProperty',[wNum oNum],'y',handles.exper.worlds{wNum}.objects{oNum}.getValue.y}};
    case 'changeSymbolicObjectLocations'
        if oNum == 0
            errordlg('Use properties to motify initial conditions.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        str = {'',''};
        
        if ischar(handles.exper.worlds{wNum}.objects{oNum}.symbolic.x)
            str{1} = handles.exper.worlds{wNum}.objects{oNum}.symbolic.x;
            str{2} = handles.exper.worlds{wNum}.objects{oNum}.symbolic.y;
            undo.param = {{'changeProperty',[wNum oNum],'x',handles.exper.worlds{wNum}.objects{oNum}.symbolic.x};
                {'changeProperty',[wNum oNum],'y',handles.exper.worlds{wNum}.objects{oNum}.symbolic.y}};
        else
            undo.param = {{'changeProperty',[wNum oNum],'x',handles.exper.worlds{wNum}.objects{oNum}.getValue.x};
                {'changeProperty',[wNum oNum],'y',handles.exper.worlds{wNum}.objects{oNum}.getValue.y}};
        end
        
        answer = inputdlg({'Symbolic expression for object X','Symbolic expression for object Y'},'Symbolic',1,str);
        if isempty(answer)
            set(guifig,'Pointer','arrow'); return
        end
        redo.param = cell(0,1);
        if ~isempty(answer{1})
            handles.exper.worlds{wNum}.objects{oNum}.x = answer{1};
            redo.param{end+1,1} = {'changeProperty',[wNum oNum],'x',answer{1}};
        end
        if ~isempty(answer{2})
            handles.exper.worlds{wNum}.objects{oNum}.y = answer{2};
            redo.param{end+1,1} = {'changeProperty',[wNum oNum],'y',answer{2}};
        end
        if size(redo.param,1)==1
            redo.param = redo.param{1};
        end
    case 'renameShape'
        if oNum == 0
            errordlg('No shape selected.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        str = get(handles.pop_shape,'string');
        val = get(handles.pop_shape,'value');
        answer = inputdlg({['New name for ' str{val}]},'New name',1,str(val));
        if isempty(answer)
            set(guifig,'Pointer','arrow'); return
        end
        if ~isvarname(answer{1})
            errordlg(['''' answer{1} ''' is an invalid variable name.'],'Error');
            set(guifig,'Pointer','arrow'); return
        end
        undo.param = {'changeProperty',[wNum oNum NaN sNum],'name',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.name};
        handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.name = answer{1};
        redo.param = {'changeProperty',[wNum oNum NaN sNum],'name',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.name};
    case 'sortShapes'
        if oNum == 0
            errordlg('No object with texture selected.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        if length(handles.exper.worlds{wNum}.objects{oNum}.texture.shapes)<3
            errordlg('Need at least two shapes other than the boundary to reorder.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        fld = cell(1,length(handles.exper.worlds{wNum}.objects{oNum}.texture.shapes)-1);
        for ndx = 1:length(fld)
            fld{ndx} = handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{ndx+1}.fullName;
        end
        ord = reorderVariables(fld,'Shapes');
        if ~all(ord==(1:length(ord)))
            [~,orig] = sort(ord);
            undo.param = {'sortObjects',wNum,oNum,orig};
            redo.param = {'sortObjects',wNum,oNum,ord};
            handles.exper.worlds{wNum}.objects{oNum}.texture.shapes(2:end) = handles.exper.worlds{wNum}.objects{oNum}.texture.shapes(ord+1);
            if sNum > 1
                handles.state.selectedShape = find(ord==sNum-1)+1;
            end
        else
            set(guifig,'Pointer','arrow'); return
        end
    case 'changeShapeProperties'
        undo.param = cell(0,1);
        redo.param = cell(0,1);
        for ndx = 1:length(inp{1})
            undo.param{end+1,1} = {'changeProperty',[wNum oNum NaN sNum],inp{1}{ndx},handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.getValue.(inp{1}{ndx})};
            handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.(inp{1}{ndx}) = inp{2}{ndx};
            redo.param{end+1,1} = {'changeProperty',[wNum oNum NaN sNum],inp{1}{ndx},handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.getValue.(inp{1}{ndx})};
        end
        if size(undo.param,1)==1
            undo.param = undo.param{1};
            redo.param = redo.param{1};
        end
    case 'changeShapeLocations'
        if ischar(handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.symbolic.x)
            undo.param = {{'changeProperty',[wNum oNum NaN sNum],'x',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.symbolic.x};
                {'changeProperty',[wNum oNum NaN sNum],'y',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.symbolic.y}};
        else
            undo.param = {{'changeProperty',[wNum oNum NaN sNum],'x',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.getValue.x};
                {'changeProperty',[wNum oNum NaN sNum],'y',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.getValue.y}};
        end
        handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.x = inp(:,1);
        handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.y = inp(:,2);
        redo.param = {{'changeProperty',[wNum oNum NaN sNum],'x',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.getValue.x};
            {'changeProperty',[wNum oNum NaN sNum],'y',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.getValue.y}};
    case 'addShapeLocation'
        if oNum == 0
            errordlg('No shape selected.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        val = get(handles.table_shapeLocations,'data');
        val = val(:,2:end);
        if ischar(handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.symbolic.x)
            undo.param = {{'changeProperty',[wNum oNum NaN sNum],'x',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.symbolic.x};
                {'changeProperty',[wNum oNum NaN sNum],'y',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.symbolic.y}};
        else
            undo.param = {{'changeProperty',[wNum oNum NaN sNum],'x',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.getValue.x};
                {'changeProperty',[wNum oNum NaN sNum],'y',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.getValue.y}};
        end
        val = [val; {'0','0'}];
        handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.x = val(:,1);
        handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.y = val(:,2);
        redo.param = {{'changeProperty',[wNum oNum NaN sNum],'x',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.getValue.x};
            {'changeProperty',[wNum oNum NaN sNum],'y',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.getValue.y}};
    case 'deleteShapeLocations'
        if oNum == 0
            errordlg('No shape selected.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        indx = get(handles.table_shapeLocations,'userdata');
        if isempty(indx)
            errordlg('No locations selected.','Error')
            set(guifig,'Pointer','arrow'); return
        end
        val = get(handles.table_shapeLocations,'data');
        val = val(:,2:end);
        if size(val,1) == 1
            errordlg('At least one location is required.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        if ischar(handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.symbolic.x)
            undo.param = {{'changeProperty',[wNum oNum NaN sNum],'x',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.symbolic.x};
                {'changeProperty',[wNum oNum NaN sNum],'y',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.symbolic.y}};
        else
            undo.param = {{'changeProperty',[wNum oNum NaN sNum],'x',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.getValue.x};
                {'changeProperty',[wNum oNum NaN sNum],'y',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.getValue.y}};
        end
        val(indx,:) = [];
        handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.x = val(:,1);
        handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.y = val(:,2);
        redo.param = {{'changeProperty',[wNum oNum NaN sNum],'x',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.getValue.x};
            {'changeProperty',[wNum oNum NaN sNum],'y',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.getValue.y}};
    case 'changeSymbolicShapeLocations'
        if oNum == 0
            errordlg('No shape selected.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        str = {'',''};
        if ischar(handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.symbolic.x)
            str{1} = handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.symbolic.x;
            str{2} = handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.symbolic.y;
            undo.param = {{'changeProperty',[wNum oNum NaN sNum],'x',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.symbolic.x};
                {'changeProperty',[wNum oNum NaN sNum],'y',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.symbolic.y}};
        else
            undo.param = {{'changeProperty',[wNum oNum NaN sNum],'x',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.getValue.x};
                {'changeProperty',[wNum oNum NaN sNum],'y',handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.getValue.y}};
        end
        answer = inputdlg({'Symbolic expression for shape X','Symbolic expression for shape Y'},'Symbolic',1,str);
        if isempty(answer)
            set(guifig,'Pointer','arrow'); return
        end
        
        redo.param = cell(0,1);
        if ~isempty(answer{1})
            handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.x = answer{1};
            redo.param{end+1,1} = {'changeProperty',[wNum oNum NaN sNum],'x',answer{1}};
        end
        if ~isempty(answer{2})
            handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.y = answer{2};
            redo.param{end+1,1} = {'changeProperty',[wNum oNum NaN sNum],'y',answer{2}};
        end
        if size(redo.param,1)==1
            redo.param = redo.param{1};
        end
    case 'changeLayout'
        layout = handles.layouts{inp};
        layout = rmfield(layout,'name');
        layout = rmfield(layout,'icon');
        handles.figs = figureLayout(handles.figs,layout);
        
        f = handles.menus.(makeVar('Layout'));
        ch = get(f,'children');
        set(ch,'checked','off');
        f(1) = handles.menus.(makeVar('Experiment layout'));
        f(2) = handles.menus.(makeVar('World layout'));
        f(3) = handles.menus.(makeVar('Texture layout'));
        f(4) = handles.menus.(makeVar('3D world view'));
        set(f,'checked','off');
        set(ch(length(ch)-inp+1),'checked','on');
        if inp < 4
            set(f(inp),'checked','on');
        end
        
        set(guifig,'Pointer','arrow'); return
    case 'builtinLayout'
        f = handles.menus.(makeVar('Layout'));
        ch = get(f,'children');
        set(ch,'checked','off');
        f = handles.menus.(makeVar('Experiment layout'));
        set(f,'checked','off');
        f = handles.menus.(makeVar('World layout'));
        set(f,'checked','off');
        f = handles.menus.(makeVar('Texture layout'));
        set(f,'checked','off');
        f = handles.menus.(makeVar('3D world view'));
        set(f,'checked','off');
        
        layout = struct;
        switch inp
            case 'experiment'
                set(ch(length(ch)),'checked','on');
                f = handles.menus.(makeVar('Experiment layout'));
                set(f,'checked','on');
                layout = handles.layouts{1};
            case 'world'
                set(ch(length(ch)-1),'checked','on');
                f = handles.menus.(makeVar('World layout'));
                set(f,'checked','on');
                layout = handles.layouts{2};
            case 'texture'
                set(ch(length(ch)-2),'checked','on');
                f = handles.menus.(makeVar('Texture layout'));
                set(f,'checked','on');
                layout = handles.layouts{3};
            case '3d'
                layout.worldDrawing = [0 0 1 1];
                f = handles.menus.(makeVar('3D world view'));
                set(f,'checked','on');
        end
        if ~strcmp(inp,'3d')
            layout = rmfield(layout,'name');
            layout = rmfield(layout,'icon');
        end
        handles.figs = figureLayout(handles.figs,layout);
        if strcmp(inp,'3d')
            rotate3d(guifig,'on');
        end
        set(guifig,'Pointer','arrow'); return
    case 'changeVariables'
        undo.param = {'changeVariables',inp{1},handles.exper.variables.(inp{1})};
        handles.exper.variables.(inp{1}) = inp{2};
        redo.param = {'changeVariables',inp{1},handles.exper.variables.(inp{1})};
    case 'saveExperiment'
        mfile = mfilename('fullpath');
        path = fileparts(mfile);
        if isempty(handles.state.fileName)
            filename = '*.mat';
        else
            filename = handles.state.fileName;
        end
        [filename, pathname] = uiputfile([path filesep '..' filesep '..' filesep 'experiments' filesep filename],'Save experiment');
        if ~ischar(filename)
            figure(guifig);
            set(guifig,'Pointer','arrow'); return
        end
        if strcmp(func2str(handles.exper.experimentCode),'undefined')
            mfileName = filename;
            if ~exist([path filesep '..' filesep '..' filesep 'experiments' filesep mfileName(1:end-4) '.m'],'file')
                handles.exper.experimentCode = str2func(mfileName(1:end-4));
                fidIn = fopen([path filesep '..' filesep '..' filesep 'defaults' filesep 'defaultVirmenCode.m']);
                fidOut = fopen([path filesep '..' filesep '..' filesep 'experiments' filesep mfileName(1:end-4) '.m'],'w'); %#ok<MCMFL>
                while 1
                    tline = fgetl(fidIn);
                    if ~ischar(tline)
                        break
                    end
                    tline = strrep(tline,'defaultVirmenCode',mfileName(1:end-4));
                    tline = strrep(tline,'%','%%');
                    fprintf(fidOut,[tline '\n']);
                end
                fclose(fidIn);
                fclose(fidOut);
            end
            edit([path filesep '..' filesep '..' filesep 'experiments' filesep mfileName(1:end-4) '.m']);
        end
        handles.exper.name = filename(1:end-4);
        exper = handles.exper; %#ok<NASGU>
        save([pathname filename],'exper');
        handles.state.fileName = filename;
        justSaved = true;
        figure(guifig);
    case 'openExperiment'
        f = handles.buttons.(makeVar('Save experiment'));
        if strcmp(get(f,'enable'),'on')
            button = questdlg('Save changes before closing?','Save','Yes','No','Cancel','Cancel');
            switch button
                case 'Yes'
                    virmenEventHandler('saveExperiment');
                    if strcmp(get(f,'enable'),'on')
                        figure(guifig);
                        set(guifig,'Pointer','arrow'); return
                    end
                case 'No'
                    % Do nothing
                case 'Cancel'
                    figure(guifig);
                    set(guifig,'Pointer','arrow'); return
            end
        end
        
        mfile = mfilename('fullpath');
        path = fileparts(mfile);
        [filename, pathname] = uigetfile([path filesep '..' filesep '..' filesep 'experiments' filesep '*.mat'],'Open experiment');
        if ~ischar(filename)
            figure(guifig);
            set(guifig,'Pointer','arrow'); return
        end
        load([pathname filename],'exper');
        handles.exper = exper; %#ok<NODEF>
        handles.exper.enableCallbacks;
        
        handles.state.fileName = filename;
        handles.state.selectedWorld = 1;
        handles.state.selectedObject = 0;
        handles.state.selectedShape = 1;
        justSaved = true;
        handles.history.position = 1;
        handles.history.states{handles.history.position}.state = handles.state;
        for ndx = 2:length(handles.history.states)
            handles.history.states{ndx}.state = [];
        end
        
        mfile = mfilename('fullpath');
        path = fileparts(mfile);
        mf = dir([path filesep '..' filesep '..' filesep 'experiments' filesep '*.m']);
        isFound = false;
        for ndx = 1:length(mf)
            f = strfind(mf(ndx).name,'.');
            if strcmp(mf(ndx).name(1:f(end)-1),func2str(handles.exper.experimentCode))
                isFound = true;
            end
        end
        if ~isFound
            fid = fopen([path filesep '..' filesep '..' filesep 'experiments' filesep func2str(handles.exper.experimentCode) '.m'],'w');
            for ndx = 1:length(handles.exper.codeText)
                fprintf(fid,'%s',handles.exper.codeText{ndx});
                fprintf(fid,'\n');
            end
            fclose(fid);
            warndlg(['Associated .m file was not found. A file ''' func2str(handles.exper.experimentCode) '.m'' was created with recovered code.'],'Warning','modal');
        end
        figure(guifig);
        
    case 'newExperiment'
        f = handles.buttons.(makeVar('Save experiment'));
        if strcmp(get(f,'enable'),'on')
            button = questdlg('Save changes before closing?','Save','Yes','No','Cancel','Cancel');
            switch button
                case 'Yes'
                    virmenEventHandler('saveExperiment');
                    if strcmp(get(f,'enable'),'on')
                        set(guifig,'Pointer','arrow'); return
                    end
                case 'No'
                    % Do nothing
                case 'Cancel'
                    set(guifig,'Pointer','arrow'); return
            end
        end
        
        handles.exper = virmenExperiment;
        handles.exper.worlds{1}.backgroundColor = handles.defaultProperties.worldBackgroundColor;
        handles.exper.worlds{1}.startLocation = handles.defaultProperties.startLocation;
        
        mfile = mfilename('fullpath');
        path = fileparts(mfile);
        fid = fopen([path filesep '..' filesep '..' filesep 'defaults' filesep 'defaultFunctions.txt']);
        txt = textscan(fid,'%s','delimiter','\t');
        txt = txt{1};
        for ndx = 1:2:length(txt)-1
            handles.exper.(txt{ndx}) = str2func(txt{ndx+1});
        end
        fclose(fid);
        
        handles.state = virmenGuiState(handles.defaultProperties);
        justSaved = true;
        handles.history.position = 1;
        handles.history.states{handles.history.position}.state = handles.state;
        for ndx = 2:length(handles.history.states)
            handles.history.states{ndx}.state = [];
        end
    case 'changeHistory'
        switch inp
            case 'undo'
                historyAct = handles.history.states{handles.history.position}.undo;
                handles.history.position = handles.history.position-1;
                handles.state = handles.history.states{handles.history.position}.state;
            case 'redo'
                historyAct = handles.history.states{handles.history.position}.redo;
                handles.history.position = handles.history.position+1;
                handles.state = handles.history.states{handles.history.position}.state;
        end
        if size(historyAct.param,1) == 1
            historyAct.param = {historyAct.param};
        end
        for j = 1:size(historyAct.param,1)
            switch historyAct.param{j}{1}
                case 'runMacro'
                    exper = historyAct.param{j}{2};
                    historyError = historyAct.param{j}{3};
                    handles.exper = exper;
                    if historyError
                        errordlg(['Cannot ' inp ' because this macro has history tracking disabled.'],'Error');
                        set(guifig,'Pointer','arrow'); return
                    end
                case 'macroSettings'
                    macroType = historyAct.param{j}{2};
                    propVal = historyAct.param{j}{3};
                    props = fieldnames(propVal);
                    for ndx = 1:length(props)
                        handles.macros.(macroType).(props{ndx}) = propVal.(props{ndx});
                    end
                case 'addShape'
                    wNum = historyAct.param{j}{2};
                    oNum = historyAct.param{j}{3};
                    sNum = historyAct.param{j}{4};
                    shape = historyAct.param{j}{5};
                    addShape(handles.exper.worlds{wNum}.objects{oNum}.texture,shape);
                    handles.exper.worlds{wNum}.objects{oNum}.texture.shapes = handles.exper.worlds{wNum}.objects{oNum}.texture.shapes([1:sNum-1 end sNum:end-1]);
                case 'deleteShape'
                    wNum = historyAct.param{j}{2};
                    oNum = historyAct.param{j}{3};
                    sNum = historyAct.param{j}{4};
                    handles.exper.worlds{wNum}.objects{oNum}.texture.shapes(sNum) = [];
                case 'addObject'
                    wNum = historyAct.param{j}{2};
                    oNum = historyAct.param{j}{3};
                    obj = historyAct.param{j}{4};
                    addObject(handles.exper.worlds{wNum},obj);
                    handles.exper.worlds{wNum}.objects = handles.exper.worlds{wNum}.objects([1:oNum-1 end oNum:end-1]);
                case 'deleteObject'
                    wNum = historyAct.param{j}{2};
                    oNum = historyAct.param{j}{3};
                    handles.exper.worlds{wNum}.objects(oNum) = [];
                case 'changeProperty'
                    indx = historyAct.param{j}{2};
                    prop = historyAct.param{j}{3};
                    val = historyAct.param{j}{4};
                    switch length(indx)
                        case 0
                            handles.exper.(prop) = val;
                        case 1
                            handles.exper.worlds{indx(1)}.(prop) = val;
                        case 2
                            handles.exper.worlds{indx(1)}.objects{indx(2)}.(prop) = val;
                        case 3
                            handles.exper.worlds{indx(1)}.objects{indx(2)}.texture.(prop) = val;
                        case 4
                            handles.exper.worlds{indx(1)}.objects{indx(2)}.texture.shapes{indx(4)}.(prop) = val;
                    end
                case 'sortObjects'
                    wNum = historyAct.param{j}{2};
                    ord = historyAct.param{j}{3};
                    handles.exper.worlds{wNum}.objects = handles.exper.worlds{wNum}.objects(ord);
                case 'sortShapes'
                    wNum = historyAct.param{j}{2};
                    oNum = historyAct.param{j}{3};
                    ord = historyAct.param{j}{4};
                    handles.exper.worlds{wNum}.objects{oNum}.texture.shapes = handles.exper.worlds{wNum}.objects{oNum}.texture.shapes(ord);
                case 'changeVariables'
                    varb = historyAct.param{j}{2};
                    val = historyAct.param{j}{3};
                    handles.exper.variables.(varb) = val;
                case 'sortWorlds'
                    ord = historyAct.param{j}{2};
                    handles.exper.worlds = handles.exper.worlds(ord);
                case 'addWorld'
                    wNum = historyAct.param{j}{2};
                    obj = historyAct.param{j}{3};
                    addWorld(handles.exper,obj);
                    handles.exper.worlds = handles.exper.worlds([1:wNum-1 end wNum:end-1]);
                case 'deleteWorld'
                    wNum = historyAct.param{j}{2};
                    handles.exper.worlds(wNum) = [];
                case 'addVariables'
                    fld = historyAct.param{j}{4};
                    for v = 1:length(historyAct.param{j}{2})
                        varb = historyAct.param{j}{2}{v};
                        val = historyAct.param{j}{3}{v};
                        handles.exper.variables.(varb) = val;
                    end
                    ord = zeros(1,length(fld));
                    currfld = fieldnames(handles.exper.variables);
                    for f = 1:length(fld)
                        ord(f) = find(cellfun(@(x)strcmp(x,fld{f}),currfld));
                    end
                    handles.exper.variables = orderfields(handles.exper.variables,ord);
                case 'deleteVariables'
                    for v = 1:length(historyAct.param{j}{2})
                        varb = historyAct.param{j}{2}{v};
                        handles.exper.variables = rmfield(handles.exper.variables,varb);
                    end
                case 'renameVariables'
                    for v = 1:length(historyAct.param{j}{2})
                        handles.exper.renameVariable(historyAct.param{j}{3}{v},['temporary_' historyAct.param{j}{2}{v}]);
                    end
                    for v = 1:length(historyAct.param{j}{2})
                        handles.exper.renameVariable(['temporary_' historyAct.param{j}{2}{v}],historyAct.param{j}{2}{v});
                    end
                case 'sortVariables'
                    ord = historyAct.param{j}{2};
                    handles.exper.variables = orderfields(handles.exper.variables,ord);
            end
        end
    case 'changeExperimentProperties'
        if ischar(inp)
            out = organizeWindows(handles.exper.windows,handles.defaultProperties);
            if isempty(out)
                set(guifig,'Pointer','arrow'); return
            end
            undo.param = {'changeProperty',[],'windows',handles.exper.windows};
            handles.exper.windows = out;
            redo.param = {'changeProperty',[],'windows',handles.exper.windows};
        else
            if isa(handles.exper.(inp{1}),'function_handle')
                undo.param = {'changeProperty',[],inp{1},handles.exper.(inp{1})};
            else
                undo.param = {'changeProperty',[],inp{1},handles.exper.getValue.(inp{1})};
            end
            handles.exper.(inp{1}) = inp{2};
            if isa(handles.exper.(inp{1}),'function_handle')
                redo.param = {'changeProperty',[],inp{1},handles.exper.(inp{1})};
            else
                redo.param = {'changeProperty',[],inp{1},handles.exper.getValue.(inp{1})};
            end
        end
    case 'selectWorld'
        if handles.state.selectedWorld ~= inp
            handles.state.selectedWorld = inp;
            handles.state.selectedObject = 0;
            handles.state.selectedShape = 1;
        end
        if strcmp(get(guifig,'selectiontype'),'open')
            virmenEventHandler('builtinLayout','world');
        end
    case 'editWorld'
        if handles.state.selectedWorld ~= inp
            handles.state.selectedWorld = inp;
            handles.state.selectedObject = 0;
            handles.state.selectedShape = 1;
        end
        guidata(guifig, handles)
        virmenEventHandler('builtinLayout','world');
    case 'renameWorld'
        if ischar(inp) && strcmp(inp,'n/a')
            inp = wNum;
        end
        answer = inputdlg({['New name for ''' handles.exper.worlds{inp}.name '''']},'New name',1,{handles.exper.worlds{inp}.name});
        if isempty(answer)
            set(guifig,'Pointer','arrow'); return
        end
        if ~isvarname(answer{1})
            errordlg(['''' answer{1} ''' is an invalid variable name.'],'Error');
            set(guifig,'Pointer','arrow'); return
        end
        undo.param = {'changeProperty',inp,'name',handles.exper.worlds{inp}.name};
        handles.exper.worlds{inp}.name = answer{1};
        redo.param = {'changeProperty',inp,'name',handles.exper.worlds{inp}.name};
    case 'sortWorlds'
        if length(handles.exper.worlds)<2
            errordlg('Need at least two worlds to reorder.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        fld = cell(1,length(handles.exper.worlds));
        for ndx = 1:length(fld)
            fld{ndx} = handles.exper.worlds{ndx}.fullName;
        end
        ord = reorderVariables(fld,'Worlds');
        if ~all(ord==(1:length(ord)))
            [~,orig] = sort(ord);
            undo.param = {'sortWorlds',orig};
            redo.param = {'sortWorlds',ord};
            handles.exper.worlds = handles.exper.worlds(ord);
            handles.state.selectedWorld = find(ord==wNum);
        else
            set(guifig,'Pointer','arrow'); return
        end
    case 'addWorld'
        world = virmenWorld;
        world.backgroundColor = handles.defaultProperties.worldBackgroundColor;
        world.transparency = handles.defaultProperties.worldTransparency;
        world.startLocation = handles.defaultProperties.startLocation;
        handles.exper = addWorld(handles.exper,world,'new');
        handles.state.selectedWorld = length(handles.exper.worlds);
        handles.state.selectedObject = 0;
        handles.state.selectedShape = 1;
        undo.param = {'deleteWorld',handles.state.selectedWorld};
        redo.param = {'addWorld',handles.state.selectedWorld,copyVirmenObject(world)};
    case 'deleteWorld'
        if length(handles.exper.worlds)==1
            errordlg('At least one world is required.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        undo.param = {'addWorld',wNum,handles.exper.worlds{wNum}};
        redo.param = {'deleteWorld',wNum};
        handles.exper.worlds(wNum) = [];
        if wNum > length(handles.exper.worlds)
            handles.state.selectedWorld = handles.state.selectedWorld - 1;
        end
        handles.state.selectedObject = 0;
        handles.state.selectedShape = 1;
    case 'addVariable'
        answer = inputdlg({'New variable name','Value for the new variable'},'New variable',1,{'',''});
        if isempty(answer)
            set(guifig,'Pointer','arrow'); return
        end
        if ~isvarname(answer{1})
            errordlg([answer{1} ' is not a valid variable name.'],'Error');
            set(guifig,'Pointer','arrow'); return
        end
        undo.param = {'deleteVariables',answer(1)};
        handles.exper.variables.(answer{1}) = answer{2};
        redo.param = {'addVariables',answer(1),answer(2),fieldnames(handles.exper.variables)};
    case 'deleteVariables'
        data = get(handles.table_variables,'data');
        if isempty(data)
            errordlg('No variables exist.','Error')
            set(guifig,'Pointer','arrow'); return
        end
        row = data(:,1);
        indx = get(handles.table_variables,'userdata');
        if isempty(indx)
            errordlg('No variables selected.','Error')
            set(guifig,'Pointer','arrow'); return
        end
        val = row(indx);
        varVals = cell(1,length(val));
        for ndx = 1:length(val)
            varVals{ndx} = handles.exper.variables.(val{ndx});
        end
        undo.param = {'addVariables',val,varVals,fieldnames(handles.exper.variables)};
        redo.param = {'deleteVariables',val};
        for ndx = 1:length(val)
            handles.exper.variables = rmfield(handles.exper.variables,val{ndx});
        end
    case 'renameVariables'
        data = get(handles.table_variables,'data');
        if isempty(data)
            errordlg('No variables exist.','Error')
            set(guifig,'Pointer','arrow'); return
        end
        row = data(:,1);
        indx = get(handles.table_variables,'userdata');
        if isempty(indx)
            errordlg('No variables selected.','Error')
            set(guifig,'Pointer','arrow'); return
        end
        val = row(indx);
        query = cell(1,length(val));
        for v = 1:length(query)
            query{v} = ['New name for ' val{v}];
        end
        answer = inputdlg(query,'Rename variables',1,val);
        if isempty(answer)
            set(guifig,'Pointer','arrow'); return
        end
        if ~isvarname(answer{1})
            errordlg([answer{1} ' is not a valid variable name.'],'Error');
            set(guifig,'Pointer','arrow'); return
        end
        undo.param = {'renameVariables',val,answer};
        redo.param = {'renameVariables',answer,val};
        for v = 1:length(query)
            handles.exper.renameVariable(val{v},['temporary_' answer{v}]);
        end
        for v = 1:length(query)
            handles.exper.renameVariable(['temporary_' answer{v}],answer{v});
        end
    case 'sortVariables'
        fld = fieldnames(handles.exper.variables);
        if length(fld)<2
            errordlg('Need at least two variables to reorder.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        ord = reorderVariables(fld,'Variables');
        if ~all(ord==(1:length(ord)))
            [~,orig] = sort(ord);
            undo.param = {'sortVariables',orig};
            redo.param = {'sortVariables',ord};
            handles.exper.variables = orderfields(handles.exper.variables,ord);
        else
            set(guifig,'Pointer','arrow'); return
        end
    case 'importWorld'
        [world, vars] = chooseObjectGui(handles.exper,'virmenWorld');
        if isempty(world)
            set(guifig,'Pointer','arrow'); return
        end
        world.enableCallbacks;
        varnames = fieldnames(vars);
        undo.param{2,1} = {'changeProperty',[],'variables',handles.exper.variables};
        for ndx = 1:length(varnames)
            handles.exper.variables.(varnames{ndx}) = vars.(varnames{ndx});
        end
        redo.param{1,1} = {'changeProperty',[],'variables',handles.exper.variables};
        handles.exper = addWorld(handles.exper,world,'copy');
        handles.state.selectedWorld = length(handles.exper.worlds);
        handles.state.selectedObject = 0;
        handles.state.selectedShape = 1;
        undo.param{1,1} = {'deleteWorld',handles.state.selectedWorld};
        redo.param{2,1} = {'addWorld',handles.state.selectedWorld,copyVirmenObject(world)};
    case 'importObject'
        [obj, vars] = chooseObjectGui(handles.exper,'virmenObject');
        if isempty(obj)
            set(guifig,'Pointer','arrow'); return
        end
        obj.enableCallbacks;
        varnames = fieldnames(vars);
        undo.param{2,1} = {'changeProperty',[],'variables',handles.exper.variables};
        for ndx = 1:length(varnames)
            handles.exper.variables.(varnames{ndx}) = vars.(varnames{ndx});
        end
        redo.param{1,1} = {'changeProperty',[],'variables',handles.exper.variables};
        handles.exper.worlds{wNum} = addObject(handles.exper.worlds{wNum},obj,'copy');
        handles.state.selectedObject = length(handles.exper.worlds{wNum}.objects);
        handles.state.selectedShape = 1;
        undo.param{1,1} = {'deleteObject',wNum,handles.state.selectedObject};
        redo.param{2,1} = {'addObject',wNum,handles.state.selectedObject,copyVirmenObject(obj)};
    case 'loadImage'
        if oNum == 0
            errordlg('No object with texture selected.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        [filename, pathname] = uigetfile('*.gif;*.jpg;*.png', 'Pick an image file');
        figure(guifig);
        if ~ischar(filename)
            set(guifig,'Pointer','arrow'); return
        end
        [img, errorString] = readTextureImageFile([pathname filename]);
        if ~isempty(errorString)
            errordlg(errorString,'Error');
            set(guifig,'Pointer','arrow'); return
        end
        if isempty(img)
            set(guifig,'Pointer','arrow'); return
        end
        undo.param = {'changeProperty',[wNum oNum],'texture',copyVirmenObject(handles.exper.worlds{wNum}.objects{oNum}.texture)};
        handles.exper.worlds{wNum}.objects{oNum}.tiling = [1 1];
        handles.exper.worlds{wNum}.objects{oNum}.texture.loadImage(img);
        if handles.state.selectedShape > length(handles.exper.worlds{wNum}.objects{oNum}.texture.shapes)
            handles.state.selectedShape = length(handles.exper.worlds{wNum}.objects{oNum}.texture.shapes);
        end
        redo.param = {'changeProperty',[wNum oNum],'texture',copyVirmenObject(handles.exper.worlds{wNum}.objects{oNum}.texture)};
    case 'importTexture'
        if oNum == 0
            errordlg('No object with texture selected.','Error');
            set(guifig,'Pointer','arrow'); return
        end
        [texture, vars] = chooseObjectGui(handles.exper,'virmenTexture');
        if isempty(texture)
            set(guifig,'Pointer','arrow'); return
        end
        texture.enableCallbacks;
        varnames = fieldnames(vars);
        undo.param{1,1} = {'changeProperty',[],'variables',handles.exper.variables};
        for ndx = 1:length(varnames)
            handles.exper.variables.(varnames{ndx}) = vars.(varnames{ndx});
        end
        redo.param{1,1} = {'changeProperty',[],'variables',handles.exper.variables};
        
        oth = zeros(0,2);
        for w = 1:length(handles.exper.worlds)
            for o = 1:length(handles.exper.worlds{w}.objects)
                if ~all([w o]==[wNum oNum]) && isequalwithequalnans(handles.exper.worlds{w}.objects{o}.texture.triangles, ...
                        handles.exper.worlds{wNum}.objects{oNum}.texture.triangles) %#ok<DISEQN>
                    oth(end+1,:) = [w o]; %#ok<AGROW>
                end
            end
        end
        undo.param{end+1,1} = {'changeProperty',[wNum oNum],'texture',copyVirmenObject(handles.exper.worlds{wNum}.objects{oNum}.texture)};
        handles.exper.worlds{wNum}.objects{oNum} = setTexture(handles.exper.worlds{wNum}.objects{oNum},texture,'copy');
        redo.param{end+1,1} = {'changeProperty',[wNum oNum],'texture',copyVirmenObject(handles.exper.worlds{wNum}.objects{oNum}.texture)};
        if ~isempty(oth)
            str = cell(1,size(oth,1));
            for ndx = 1:size(oth,1)
                str{ndx} = handles.exper.worlds{oth(ndx,1)}.objects{oth(ndx,2)}.fullName;
            end
            [indx,ok] = listdlg('ListString',str,'InitialValue',1:length(str),'ListSize',[250 150], ...
                'Name','Select objects','PromptString','Change other objects with identical texture:');
            if ok > 0
                for ndx = 1:length(indx)
                    undo.param{end+1,1} = {'changeProperty',oth(indx(ndx),:),'texture',copyVirmenObject(handles.exper.worlds{oth(indx(ndx),1)}.objects{oth(indx(ndx),2)}.texture)};
                    setTexture(handles.exper.worlds{oth(indx(ndx),1)}.objects{oth(indx(ndx),2)},handles.exper.worlds{wNum}.objects{oNum}.texture,'copy');
                    redo.param{end+1,1} = {'changeProperty',oth(indx(ndx),:),'texture',copyVirmenObject(handles.exper.worlds{oth(indx(ndx),1)}.objects{oth(indx(ndx),2)}.texture)};
                end
            end
        end
        undo.param = undo.param([2:end 1]);
        handles.state.selectedShape = 1;
    case 'run'
        figs = findall(0,'type','figure','visible','on');
        set(figs,'visible','off');
        if strcmp(func2str(handles.exper.experimentCode),'undefined')
            handles.exper.experimentCode = @defaultVirmenCode;
            isDefault = true;
        else
            isDefault = false;
        end
        drawnow
        err = virmenEngine(handles.exper);
        if isDefault
            handles.exper.experimentCode = @undefined;
        end
        set(figs,'visible','on');
        figure(guifig);
        set(guifig,'Pointer','arrow');
        if isstruct(err)
            errordlg(['User function ' err.stack(1).name ' generated an error on line ' num2str(err.stack(1).line) ': ' err.message],'Error');
            if length(err.stack) > 1
                err.stack(2:end) = [];
            end
            error(err);
        end
        return
    case 'closeProgram'
        f = handles.buttons.(makeVar('Save experiment'));
        if strcmp(get(f,'enable'),'on')
            button = questdlg('Save changes before closing?','Save','Yes','No','Cancel','Cancel');
            switch button
                case 'Yes'
                    virmenEventHandler('saveExperiment');
                    if strcmp(get(f,'enable'),'on')
                        set(guifig,'Pointer','arrow'); return
                    end
                case 'No'
                    % Do nothing
                case 'Cancel'
                    set(guifig,'Pointer','arrow'); return
            end
        end
        delete(handles.mainFigure);
        varfig = findall(0,'name','ViRMEn variables');
        delete(varfig);
        aboutfig = findall(0,'name','About ViRMEn');
        delete(aboutfig);
        return
    case 'listVariables'
        variablesGui;
    case 'export'
        switch inp
            case 'World 2D'
                world = handles.exper.worlds{wNum};
                figure('Name',world.fullName);
                world.draw2D;
                view(gca,2)
                xl = handles.state.worldXLim;
                yl = handles.state.worldYLim;
                set(gca,'units','pixels');
                pos = get(gca,'position');
                set(gca,'units','normalized');
                ar = pos(4)/pos(3);
                if range(yl)/range(xl) > ar
                    xl = mean(xl) + (xl-mean(xl))/range(xl)*range(yl)/ar;
                else
                    yl = mean(yl) + (yl-mean(yl))/range(yl)*range(xl)*ar;
                end
                set(gca,'xlim',xl);
                set(gca,'ylim',yl);
                axis(gca,'equal')
            case 'World 3D'
                world = handles.exper.worlds{wNum};
                figure('Name',world.fullName);
                world.draw3D;
                axis(gca,'equal')
                axis(gca,'tight')
                set(gca,'color',world.backgroundColor);
            case 'World wireframe'
                world = handles.exper.worlds{wNum};
                figure('Name',world.fullName);
                [~, he, hp] = world.draw2D;
                delete([he hp]);
                view(gca,3)
                axis(gca,'equal')
                axis(gca,'tight')
            case 'Texture'
                if oNum == 0
                    errordlg('No object with texture selected.','Error');
                    set(guifig,'Pointer','arrow'); return
                end
                
                texture = handles.exper.worlds{wNum}.objects{oNum}.texture;
                
                figure('name',texture.fullName);
                h = texture.draw;
                
                if handles.state.showTriangulation == 0
                    set(h,'edgecolor','none');
                else
                    set(h,'edgecolor',handles.state.triangulationColor);
                end
                
                view(gca,2);
                
                if isempty(handles.state.textureXLim)
                    xl = [-0.1*texture.width 1.1*texture.width];
                    yl = [-0.1*texture.height 1.1*texture.height];
                else
                    xl = handles.state.textureXLim;
                    yl = handles.state.textureYLim;
                end
                set(gca,'xlim',xl);
                set(gca,'ylim',yl);
                set(gca,'units','pixels');
                pos = get(gca,'position');
                set(gca,'units','normalized');
                ar = pos(4)/pos(3);
                if range(yl)/range(xl) > ar
                    xl = mean(xl) + (xl-mean(xl))/range(xl)*range(yl)/ar;
                else
                    yl = mean(yl) + (yl-mean(yl))/range(yl)*range(xl)*ar;
                end
                set(gca,'xlim',xl);
                set(gca,'ylim',yl);
                set(gca,'color',handles.exper.worlds{wNum}.backgroundColor);
                set(gca,'box','on');
        end
    case 'manual'
        mfile = mfilename('fullpath');
        path = fileparts(mfile);
        open([path filesep '..' filesep 'documentation' filesep 'ViRMEn manual.pdf'])
        set(guifig,'Pointer','arrow'); return
    case 'about'
        virmenAboutWindow;
        set(guifig,'Pointer','arrow'); return
    case 'resizeFigure'
        %
end

f = handles.buttons.(makeVar('Save experiment'));

if handles.historyBool.(type) && ~(isempty(undo.param) && isequal(oldState,handles.state))
    undo.string = label2string(type);
    redo.string = undo.string;
    
    if handles.history.position == length(handles.history.states)
        handles.history.states = handles.history.states([2:end 1]);
        handles.history.states{end}.state = handles.state;
        handles.history.states{end}.undo = undo;
        handles.history.states{end-1}.redo = redo;
    else
        handles.history.position = handles.history.position + 1;
        handles.history.states{handles.history.position}.state = handles.state;
        handles.history.states{handles.history.position}.undo = undo;
        if handles.history.position > 1
            handles.history.states{handles.history.position-1}.redo = redo;
        end
        for ndx = handles.history.position+1:length(handles.history.states)
            handles.history.states{ndx}.state = [];
        end
    end
    set(f,'enable','on');
end

if handles.history.position > 0
    handles.history.states{handles.history.position}.state = handles.state;
end

if strcmp(type,'changeHistory')
    set(f,'enable','on');
end

if justSaved
    set(f,'enable','off');
end

f = handles.buttons.(makeVar('Undo'));
if handles.history.position == 1
    set(f,'enable','off');
    set(f,'tooltipstring','Undo');
else
    set(f,'enable','on');
    set(f,'tooltipstring',['Undo ' handles.history.states{handles.history.position}.undo.string]);
end
g = handles.menus.(makeVar('Undo'));
set(g,'label',get(f,'tooltipstring'),'userdata',get(f,'tooltipstring'));

f = handles.buttons.(makeVar('Redo'));
if handles.history.position == length(handles.history.states) || isempty(handles.history.states{handles.history.position+1}.state)
    set(f,'enable','off');
    set(f,'tooltipstring','Redo');
else
    set(f,'enable','on');
    set(f,'tooltipstring',['Redo ' handles.history.states{handles.history.position}.redo.string]);
end
g = handles.menus.(makeVar('Redo'));
set(g,'label',get(f,'tooltipstring'),'userdata',get(f,'tooltipstring'));

guidata(guifig, handles);
updateFigures(type);

if ~strcmp(type,'changeVariables') && ~strcmp(type,'listVariables')
    set(0,'currentfigure',guifig);
end
set(guifig,'Pointer','arrow');

function str = label2string(type)

str = type;
for ndx = 65:90
    str = strrep(str,char(ndx),[' ' char(ndx+32)]);
end
f = strfind(str,' ');

if ~isempty(f)
    f = f(1);
    if strcmp(str(f-1),'e')
        str = [str(1:f-2) 'ing' str(f:end)];
    elseif strcmp(str(f-1),'n')
        str = [str(1:f-1) 'ning' str(f:end)];
    else
        str = [str(1:f-1) 'ing' str(f:end)];
    end
end

function str = makeVar(str)

str = str(regexp(str,'[A-Za-z]'));

function ptr = zoomPointer

ptr = [   NaN   NaN   NaN   NaN     1     1     1     1   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
    NaN   NaN     1     1   NaN     2   NaN     2     1     1   NaN   NaN   NaN   NaN   NaN   NaN
    NaN     1     2   NaN     2     1     1   NaN     2   NaN     1   NaN   NaN   NaN   NaN   NaN
    NaN     1   NaN     2   NaN     1     1     2   NaN     2     1   NaN   NaN   NaN   NaN   NaN
    1   NaN     2   NaN     2     1     1   NaN     2   NaN     2     1   NaN   NaN   NaN   NaN
    1     2     1     1     1     1     1     1     1     1   NaN     1   NaN   NaN   NaN   NaN
    1   NaN     1     1     1     1     1     1     1     1     2     1   NaN   NaN   NaN   NaN
    1     2   NaN     2   NaN     1     1     2   NaN     2   NaN     1   NaN   NaN   NaN   NaN
    NaN     1     2   NaN     2     1     1   NaN     2   NaN     1   NaN   NaN   NaN   NaN   NaN
    NaN     1   NaN     2   NaN     1     1     2   NaN     2     1     2   NaN   NaN   NaN   NaN
    NaN   NaN     1     1     2   NaN     2   NaN     1     1     1     1     2   NaN   NaN   NaN
    NaN   NaN   NaN   NaN     1     1     1     1   NaN     2     1     1     1     2   NaN   NaN
    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     2     1     1     1     2   NaN
    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     2     1     1     1     2
    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     2     1     1     1
    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     2     1     2];