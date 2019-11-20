function updateFigures(type)

global guifig;
handles = guidata(guifig);

if ~isfield(handles,'figs') % figures not yet created
    return
end

handles.exper.updateNames;
set(0,'currentfigure',guifig);

figs = fieldnames(handles.figNames);
if strcmp(type,'LAYOUT')
    for ndx = 1:length(figs)
        if strcmp(get(handles.figs.(figs{ndx}),'visible'),'on')
            handles.bools(handles.evtNames.(type),handles.figNames.(figs{ndx})) = 1-handles.changedFigs(handles.figNames.(figs{ndx}));
        else
            handles.bools(handles.evtNames.(type),handles.figNames.(figs{ndx})) = 0;
        end
    end
end

for ndx = 1:length(figs)
    if strcmp(get(handles.figs.(figs{ndx}),'visible'),'on') % figure visible
        handles.changedFigs(handles.figNames.(figs{ndx})) = 1;
    elseif handles.bools(handles.evtNames.(type),handles.figNames.(figs{ndx})) == 1 % figure not visible, needs to be plotted
        handles.changedFigs(handles.figNames.(figs{ndx})) = 0;
    end
end

axisBox = [0 0 1 1];

wNum = handles.state.selectedWorld;
oNum = handles.state.selectedObject;
sNum = handles.state.selectedShape;

% TEXTURESKETCH
if strcmp(get(handles.figs.textureSketch,'visible'),'on') && handles.bools(handles.evtNames.(type),handles.figNames.textureSketch)
    handles.changedFigs(handles.figNames.textureSketch) = 1;
    
    ax = findobj(handles.figs.textureSketch,'type','axes');
    if ~isempty(ax)
        set(guifig,'currentaxes',findobj(handles.figs.textureSketch,'type','axes'));
        cla(ax)
    else
        ax = axes('parent',handles.figs.textureSketch,'outerposition',axisBox);
    end
    
    if oNum == 0 % Starting position, rather that object selected
        set(handles.figs.textureSketch,'title','Texture sketch');
        text(0,0,'No object with texture selected','horizontalalignment','center','fontsize',14,'parent',ax);
        set(ax,'xlim',[-1 1]);
        set(ax,'ylim',[-1 1]);
        axis(ax,'off');
    else
        texture = handles.exper.worlds{wNum}.objects{oNum}.texture;
        set(handles.figs.textureSketch,'title',['Sketch of ' handles.exper.worlds{wNum}.objects{oNum}.texture.indexedName ' for ' handles.exper.worlds{wNum}.objects{oNum}.indexedName]);
        axis(ax,'on');
        set(guifig,'handlevisibility','on');
        h = texture.sketch;
        set(guifig,'handlevisibility','off');
        set(h(sNum),'markersize',10,'color','r');
        set(ax,'children',[h(sNum) h([1:sNum-1 sNum+1:end])]);
        for ndx = 1:length(h)
            set(h(ndx),'buttondownfcn',['virmenEventHandler(''selectShape'',{' num2str(ndx) ',''' get(h(ndx),'marker') ''',' num2str(get(h(ndx),'markersize')) '});']);
        end
        
        if isempty(handles.state.textureXLim)
            handles.state.textureXLim = [-0.025*texture.width 1.025*texture.width];
            handles.state.textureYLim = [-0.025*texture.height 1.025*texture.height];
        end
        set(ax,'xlim',handles.state.textureXLim);
        set(ax,'ylim',handles.state.textureYLim);
        
        xl = handles.state.textureXLim;
        yl = handles.state.textureYLim;
        set(ax,'units','pixels');
        pos = get(ax,'position');
        set(ax,'units','normalized');
        ar = pos(4)/pos(3);
        if range(yl)/range(xl) > ar
            xl = mean(xl) + (xl-mean(xl))/range(xl)*range(yl)/ar;
        else
            yl = mean(yl) + (yl-mean(yl))/range(yl)*range(xl)*ar;
        end
        set(ax,'xlim',xl);
        set(ax,'ylim',yl);
        
        set(ax,'buttondownfcn','virmenEventHandler(''zoomOnTexture'',[]);');
    end
end

% TEXTUREDRAWING
if strcmp(get(handles.figs.textureDrawing,'visible'),'on') && handles.bools(handles.evtNames.(type),handles.figNames.textureDrawing)
    ax = findobj(handles.figs.textureDrawing,'type','axes');
    if ~isempty(ax)
        set(guifig,'currentaxes',findobj(handles.figs.textureDrawing,'type','axes'));
        cla(ax)
    else
        ax = axes('parent',handles.figs.textureDrawing,'outerposition',axisBox);
    end
    
    if oNum == 0
        set(handles.figs.textureDrawing,'title','Texture drawing');
        text(0,0,'No object with texture selected','horizontalalignment','center','fontsize',14,'color',1-handles.exper.worlds{wNum}.backgroundColor,'parent',ax);
        set(ax,'xlim',[-1 1]);
        set(ax,'ylim',[-1 1]);
    else
        texture = handles.exper.worlds{wNum}.objects{oNum}.texture;
        set(guifig,'handlevisibility','on');
        h = texture.draw;
        hold on
        set(guifig,'handlevisibility','off');
        
        if handles.state.showTriangulation == 0
            set(h,'edgecolor','none');
        else
            set(h,'edgecolor',handles.state.triangulationColor);
        end
        
        set([ax h],'buttondownfcn','virmenEventHandler(''zoomOnTexture'',[]);');
        num = length(find(~isnan(texture.triangles.cdata(texture.triangles.triangulation(:,1),1))));
        set(handles.figs.textureDrawing,'title',['Drawing of ' handles.exper.worlds{wNum}.objects{oNum}.texture.indexedName ' for ' handles.exper.worlds{wNum}.objects{oNum}.indexedName ': ' num2str(num) ' visible triangles']);
        
        view(ax,2);
        
        if isempty(handles.state.textureXLim)
            handles.state.textureXLim = [-0.025*texture.width 1.025*texture.width];
            handles.state.textureYLim = [-0.025*texture.height 1.025*texture.height];
        end
        set(ax,'xlim',handles.state.textureXLim);
        set(ax,'ylim',handles.state.textureYLim);
        
        xl = handles.state.textureXLim;
        yl = handles.state.textureYLim;
        set(ax,'units','pixels');
        pos = get(ax,'position');
        set(ax,'units','normalized');
        ar = pos(4)/pos(3);
        if range(yl)/range(xl) > ar
            xl = mean(xl) + (xl-mean(xl))/range(xl)*range(yl)/ar;
        else
            yl = mean(yl) + (yl-mean(yl))/range(yl)*range(xl)*ar;
        end
        set(ax,'xlim',xl);
        set(ax,'ylim',yl);
    end
    
    set(ax,'xtick',[],'ytick',[],'color',handles.exper.worlds{wNum}.backgroundColor);
    set(ax,'box','on');
end

% WORLDSKETCH
if strcmp(get(handles.figs.worldSketch,'visible'),'on') && handles.bools(handles.evtNames.(type),handles.figNames.worldSketch)
    ax = findobj(handles.figs.worldSketch,'type','axes');
    if ~isempty(ax)
        ax = findobj(handles.figs.worldSketch,'type','axes');
        set(guifig,'currentaxes',ax);
    else
        ax = axes('parent',handles.figs.worldSketch,'outerposition',axisBox);
    end
    world = handles.exper.worlds{wNum};
    
    if strcmp(type,'selectObject')
        ud = get(ax,'userdata');
        h = ud{1};
        he = ud{2};
        hp = ud{3};
    else
        cla(ax)
        set(guifig,'handlevisibility','on');
        [h, he, hp] = world.draw2D;
        set(guifig,'handlevisibility','off');
        set(ax,'userdata',{h, he, hp});
        view(ax,2);
        set(handles.figs.worldSketch,'title',[world.indexedName ' sketch'])
        for ndx = 1:length(h)
            set(h(ndx),'zdata',zeros(size(get(h(ndx),'zdata'))));
            set(he(ndx),'zdata',zeros(size(get(he(ndx),'zdata'))));
        end
        set(hp,'zdata',zeros(size(get(hp,'zdata'))));
        for ndx = 1:length(h)
            set(h(ndx),'buttondownfcn',['virmenEventHandler(''selectObject'',{' num2str(ndx) ',''' get(h(ndx),'marker') ''',' num2str(get(h(ndx),'markersize')) '});']);
            set(he(ndx),'buttondownfcn',['virmenEventHandler(''selectObject'',{' num2str(ndx) ',''' get(he(ndx),'marker') ''',' num2str(get(he(ndx),'markersize')) '});']);
        end
        set(hp,'buttondownfcn','virmenEventHandler(''selectObject'',{0,''none'',1});');
    end
    
    set([h he hp],'color','k');
    if oNum == 0
        set(hp,'color','r');
        set(ax,'children',[hp h he]);
    else
        set(h(oNum),'color','r');
        set(he(oNum),'color','r');
        set(ax,'children',[h(oNum) he(oNum) h([1:oNum-1 oNum+1:end]) he([1:oNum-1 oNum+1:end]) hp]);
    end
    
    
    xl = handles.state.worldXLim;
    yl = handles.state.worldYLim;
    set(ax,'units','pixels');
    pos = get(ax,'position');
    set(ax,'units','normalized');
    ar = pos(4)/pos(3);
    if range(yl)/range(xl) > ar
        xl = mean(xl) + (xl-mean(xl))/range(xl)*range(yl)/ar;
    else
        yl = mean(yl) + (yl-mean(yl))/range(yl)*range(xl)*ar;
    end
    set(ax,'xlim',xl);
    set(ax,'ylim',yl);
    
    set(ax,'buttondownfcn','virmenEventHandler(''zoomOnWorld'',[]);');
end

% WORLDDRAWING
if (strcmp(get(handles.figs.worldDrawing,'visible'),'on') && handles.bools(handles.evtNames.(type),handles.figNames.worldDrawing)) && handles.state.showWireframe==0
    ax = findobj(handles.figs.worldDrawing,'type','axes');
    if ~isempty(ax)
        set(guifig,'currentaxes',findobj(handles.figs.worldDrawing,'type','axes'));
    else
        ax = axes('parent',handles.figs.worldDrawing,'outerposition',axisBox);
        set(ax,'view',[-45 45]);
    end
    
    world = handles.exper.worlds{wNum};
    cla(ax)
    
    v = get(ax,'view');
    set(guifig,'handlevisibility','on');
    h = world.draw3D;
    set(guifig,'handlevisibility','off');
    set(ax,'view',v);
    set(handles.figs.worldDrawing,'title',[world.indexedName ' drawing: ' num2str(size(get(h,'faces'),1)) ' triangles, ' ...
        num2str(size(get(h,'vertices'),1)) ' vertices']);
    axis(ax,'equal');
    axis(ax,'tight');
    axis(ax,'off');
    set(ax,'color',world.backgroundColor);
    set(handles.figs.worldDrawing,'backgroundcolor',world.backgroundColor);
    set(handles.figs.worldDrawing,'foregroundcolor',1-world.backgroundColor);
    set(h,'buttondownfcn','virmenEventHandler(''builtinLayout'',''3d'')');
end
if (strcmp(get(handles.figs.worldDrawing,'visible'),'on') && ...
        (handles.bools(handles.evtNames.(type),handles.figNames.worldSketch)) || handles.bools(handles.evtNames.(type),handles.figNames.worldDrawing)) ...
        && handles.state.showWireframe == 1
    ax = findobj(handles.figs.worldDrawing,'type','axes');
    if ~isempty(ax)
        set(guifig,'currentaxes',findobj(handles.figs.worldDrawing,'type','axes'));
    else
        ax = axes('parent',handles.figs.worldDrawing,'outerposition',axisBox);
        set(ax,'view',[-45 45]);
    end
    
    world = handles.exper.worlds{wNum};
    
    if strcmp(type,'selectObject')
        h = get(ax,'userdata');
    else
        cla(ax)
        
        set(handles.figs.worldDrawing,'title',[world.indexedName ' drawing'])
        
        v = get(ax,'view');
        set(guifig,'handlevisibility','on');
        [h, he, hp] = world.draw2D;
        set(guifig,'handlevisibility','off');
        delete([he hp]);
        set(ax,'view',v,'userdata',h);
        
        for ndx = 1:length(h)
            set(h(ndx),'buttondownfcn',['virmenEventHandler(''selectObject'',{' num2str(ndx) ',''' get(h(ndx),'marker') ''',' num2str(get(h(ndx),'markersize')) '});']);
        end
    end
    
    set(h,'color','k');
    if oNum > 0
        set(h(oNum),'color','r');
        set(ax,'children',[h(oNum) h([1:oNum-1 oNum+1:end])]);
    end
    
    axis(ax,'equal');
    axis(ax,'tight');
    axis(ax,'off');
    set(ax,'color','w');
    set(handles.figs.worldDrawing,'backgroundcolor','w');
    set(handles.figs.worldDrawing,'foregroundcolor','k');
end

% OBJECTPROPERTIESGUI
if strcmp(get(handles.figs.objectProperties,'visible'),'on') && handles.bools(handles.evtNames.(type),handles.figNames.objectProperties)
    set(guifig,'currentaxes',handles.axes_texture);
    cla(handles.axes_texture)
    
    str = {'Initial conditions'};
    for ndx = 1:length(handles.exper.worlds{wNum}.objects)
        str{end+1} = handles.exper.worlds{wNum}.objects{ndx}.indexedName; %#ok<AGROW>
    end
    set(handles.pop_object,'string',str,'value',oNum+1);
    
    set(handles.table_objectProperties,'columnname',{'Property','Value'},'rowname',{},'columneditable',[false true]);
    set(handles.table_objectLocations,'enable','on');
    set(handles.table_objectLocations,'columnname',{'#','X','Y'},'rowname',{},'columneditable',[false true true]);
    if oNum == 0
        tt = get(handles.axes_texture,'title');
        set(tt,'string','');
        set(handles.axes_texture,'xlim',[-1 1]);
        set(handles.axes_texture,'ylim',[-1 1]);
        axis(handles.axes_texture,'off');
        set(handles.edit_tilingVertical,'enable','off','string','n/a');
        set(handles.edit_tilingHorizontal,'enable','off','string','n/a');
        set(handles.edit_edgeRadius,'enable','off','string','n/a');
        data = {'X','Y','Z','rotation','backgroundR','backgroundG','backgroundB','enableTransparency'}';
        data(:,2) = [handles.exper.worlds{wNum}.getValue.startLocation handles.exper.worlds{wNum}.getValue.backgroundColor handles.exper.worlds{wNum}.getValue.transparency]';
        set(handles.table_objectProperties,'data',data);
        set(handles.table_objectLocations,'enable','off','data',{});
        
    else
        obj = handles.exper.worlds{wNum}.objects{oNum};
        set(guifig,'handlevisibility','on');
        hand = obj.texture.draw;
        set(guifig,'handlevisibility','off');
        set(hand,'buttondownfcn','virmenEventHandler(''builtinLayout'',''texture'')');
        set(handles.axes_texture,'buttondownfcn','virmenEventHandler(''builtinLayout'',''texture'')');
        set(hand,'edgecolor','none');
        view(handles.axes_texture,2)
        axis(handles.axes_texture,'equal');
        set(handles.axes_texture,'xlim',[-obj.texture.width*.05 obj.texture.width*1.05]);
        set(handles.axes_texture,'ylim',[-obj.texture.height*.05 obj.texture.height*1.05]);
        set(handles.axes_texture,'xtick',[],'ytick',[],'box','on', 'color', handles.exper.worlds{wNum}.backgroundColor, ...
            'xcolor',handles.exper.worlds{wNum}.backgroundColor,'ycolor',handles.exper.worlds{wNum}.backgroundColor);
        tt = get(handles.axes_texture,'title');
        set(tt,'string',obj.texture.indexedName,'fontsize',10);
        set(tt,'buttondownfcn',['virmenEventHandler(''renameObject'',' num2str(oNum) ');'],'interpreter','none')
        set(handles.edit_tilingVertical,'enable','on','string',obj.getValue.tiling{1});
        set(handles.edit_tilingHorizontal,'enable','on','string',obj.getValue.tiling{2});
        set(handles.edit_edgeRadius,'enable','on','string',obj.getValue.edgeRadius);
        
        vprop = properties('virmenObject');
        props = properties(handles.exper.worlds{wNum}.objects{oNum});
        for ndx = length(props):-1:1
            if ~isempty(find(cellfun(@(x)strcmp(x,props{ndx}),vprop),1))
                props(ndx) = [];
            end
        end
        
        data = props;
        for ndx = 1:length(props)
            data(ndx,2) = handles.exper.worlds{wNum}.objects{oNum}.getValue.(props{ndx});
        end
        if isempty(data)
            data = {};
        end
        set(handles.table_objectProperties,'data',data);
        xStr = handles.exper.worlds{wNum}.objects{oNum}.getValue.x;
        yStr = handles.exper.worlds{wNum}.objects{oNum}.getValue.y;
        if length(xStr)==1 && length(yStr)>1
            xStr = repmat(xStr,length(yStr),1);
        elseif length(xStr)>1 && length(yStr)==1
            yStr = repmat(yStr,length(xStr),1);
        end
        set(handles.table_objectLocations,'data',[cellfun(@num2str,num2cell(1:numel(xStr)),'uniformoutput',false)' xStr(:) yStr(:)]);
    end
end

% SHAPEPROPERTIESGUI
if strcmp(get(handles.figs.shapeProperties,'visible'),'on') && handles.bools(handles.evtNames.(type),handles.figNames.shapeProperties)
    set(handles.table_shapeProperties,'columnname',{'Property','Value'},'rowname',{},'columneditable',[false true]);
    set(handles.table_shapeProperties,'enable','on');
    set(handles.table_shapeLocations,'columnname',{'#','X','Y'},'rowname',{},'columneditable',[false true true]);
    if oNum == 0
        set(handles.pop_shape,'value',1,'string','No object with texture selected','enable','off');
    else
        str = {};
        for ndx = 1:length(handles.exper.worlds{wNum}.objects{oNum}.texture.shapes)
            str{end+1} = handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{ndx}.indexedName; %#ok<AGROW>
        end
        set(handles.pop_shape,'string',str,'value',sNum,'enable','on');
    end
    if oNum == 0 || sNum == 1
        set(handles.table_shapeProperties,'Data',{},'enable','off');
        set(handles.table_shapeLocations,'Data',{},'enable','off');
    else
        set(handles.table_shapeProperties,'enable','on');
        set(handles.table_shapeLocations,'enable','on');
        
        vprop = properties('virmenShape');
        props = properties(handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum});
        for ndx = length(props):-1:1
            if ~isempty(find(cellfun(@(x)strcmp(x,props{ndx}),vprop),1))
                props(ndx) = [];
            end
        end
        
        data = props;
        for ndx = 1:length(props)
            data(ndx,2) = handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.getValue.(props{ndx});
        end
        if isempty(data)
            data = {};
        end
        set(handles.table_shapeProperties,'data',data);
        xStr = handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.getValue.x;
        yStr = handles.exper.worlds{wNum}.objects{oNum}.texture.shapes{sNum}.getValue.y;
        if length(xStr)==1 && length(yStr)>1
            xStr = repmat(xStr,length(yStr),1);
        elseif length(xStr)>1 && length(yStr)==1
            yStr = repmat(yStr,length(xStr),1);
        end
        
        set(handles.table_shapeLocations,'data',[cellfun(@num2str,num2cell(1:numel(xStr)),'uniformoutput',false)' xStr(:) yStr(:)]);
    end
end

% VARIABLESTABLE
h = {};
if strcmp(get(handles.figs.variablesTable,'visible'),'on') && handles.bools(handles.evtNames.(type),handles.figNames.variablesTable)
    h{end+1} = handles;
end
varfig = findall(0,'name','ViRMEn variables');
if ~isempty(varfig)
    h{end+1} = guidata(varfig);
end
for ndx = 1:length(h)
    set(h{ndx}.table_variables,'columnname',{'Variable','Value'},'rowname',{}, ...
        'columneditable',[false true]);
    props = fieldnames(handles.exper.variables);
    data = props;
    for i = 1:length(data)
        data{i,2} = handles.exper.variables.(props{i});
    end
    if isempty(data)
        data = {};
    end
    set(h{ndx}.table_variables,'data',data);
end

% EXPERIMENTPROPERTIES
if strcmp(get(handles.figs.experimentProperties,'visible'),'on') && handles.bools(handles.evtNames.(type),handles.figNames.experimentProperties)
    mfile = mfilename('fullpath');
    path = fileparts(mfile);
    
    % Experiment name
    set(handles.text_experimentName,'string',handles.exper.name);
    
    % Movement
    mf = dir([path filesep '..' filesep '..' filesep 'movements' filesep '*.m']);
    switch computer
        case 'PCWIN'
            mf = [mf; dir([path filesep '..' filesep '..' filesep 'movements' filesep '*.mexw32'])];
        case 'PCWIN64'
            mf = [mf; dir([path filesep '..' filesep '..' filesep 'movements' filesep '*.mexw64'])];
        case 'MACI64'
            mf = [mf; dir([path filesep '..' filesep '..' filesep 'movements' filesep '*.mexmaci64'])];
    end
    str = cell(1,length(mf));
    val = 0;
    for ndx = 1:length(mf)
        f = strfind(mf(ndx).name,'.');
        str{ndx} = mf(ndx).name(1:f(end)-1);
        if strcmp(mf(ndx).name(1:f(end)-1),func2str(handles.exper.movementFunction))
            val = ndx;
        end
    end
    if val == 0
        val = 1;
        handles.exper.movementFunction = str2func(str{1});
    end
    set(handles.pop_movements,'string',str,'value',val);
    
    % Transformation
    mf = dir([path filesep '..' filesep '..' filesep 'transformations' filesep '*.m']);
    switch computer
        case 'PCWIN'
            mf = [mf; dir([path filesep '..' filesep '..' filesep 'transformations' filesep '*.mexw32'])];
        case 'PCWIN64'
            mf = [mf; dir([path filesep '..' filesep '..' filesep 'transformations' filesep '*.mexw64'])];
        case 'MACI64'
            mf = [mf; dir([path filesep '..' filesep '..' filesep 'transformations' filesep '*.mexmaci64'])];
    end
    str = cell(1,length(mf));
    val = 0;
    for ndx = 1:length(mf)
        f = strfind(mf(ndx).name,'.');
        str{ndx} = mf(ndx).name(1:f(end)-1);
        if strcmp(mf(ndx).name(1:f(end)-1),func2str(handles.exper.transformationFunction))
            val = ndx;
        end
    end
    if val == 0
        val = 1;
        handles.exper.transformationFunction = str2func(str{1});
    end
    set(handles.pop_transformations,'string',str,'value',val);
    
    % Experiment code
    mf = dir([path filesep '..' filesep '..' filesep 'experiments' filesep '*.m']);
    str = cell(1,length(mf));
    val = 0;
    for ndx = 1:length(mf)
        f = strfind(mf(ndx).name,'.');
        str{ndx} = mf(ndx).name(1:f(end)-1);
        if strcmp(mf(ndx).name(1:f(end)-1),func2str(handles.exper.experimentCode))
            val = ndx;
        end
    end
    if val == 0
        val = 1;
        str = [{'[Default code]'} str];
        handles.exper.experimentCode = @undefined;
    end
    set(handles.pop_experimentCode,'string',str,'value',val);
end


% WORLDSMENU
if strcmp(get(handles.figs.worldsMenu,'visible'),'on') && handles.bools(handles.evtNames.(type),handles.figNames.worldsMenu)
    delete(findobj(handles.figs.worldsMenu,'type','axes'));
    currUnits = get(handles.figs.worldsMenu,'units');
    set(handles.figs.worldsMenu,'units','pixels');
    pos = get(handles.figs.worldsMenu,'position');
    set(handles.figs.worldsMenu,'units',currUnits);
    
    [hg, wd] = virmenArrangeGrid(pos(4),pos(3),length(handles.exper.worlds));
    
    scale = .95;
    
    ndx = 1;
    xs = [];
    ys = [];
    zs = [];
    axs = [];
    for y = hg:-1:1
        for x = 1:wd
            ax = axes('parent',handles.figs.worldsMenu,'outerposition',[(x-1)/wd+(1-scale)/(2*wd) (y-1)/hg+(1-scale)/(2*hg) scale/wd scale/hg]);
            if ndx <= length(handles.exper.worlds)
                set(guifig,'handlevisibility','on');
                [h, he, hp] = handles.exper.worlds{ndx}.draw2D;
                set(guifig,'handlevisibility','off');
                view(ax,3);
                if ~isempty(h)
                    delete([he hp]);
                    he = [];
                    hp = [];
                end
                
                axis(ax,'tight');
                axis(ax,'equal');
                
                if ndx == wNum
                    set(ax,'xcolor','r','ycolor','r','zcolor','r');
                end
                set(ax,'fontsize',10,'fontweight','bold');
                tt = get(ax,'title');
                set(tt,'string',handles.exper.worlds{ndx}.indexedName);
                xl = get(ax,'ylabel');
                set(xl,'string','Edit');
                set(xl,'backgroundcolor',[0.552941176470588   0.694117647058824   0.674509803921569]);
                set(ax,'xtick',[],'ytick',[],'ztick',[],'box','off');
                set(ax,'buttondownfcn',['virmenEventHandler(''selectWorld'',' num2str(ndx) ');'])
                set([h he hp],'buttondownfcn',['virmenEventHandler(''selectWorld'',' num2str(ndx) ');'])
                set(tt,'buttondownfcn',['virmenEventHandler(''renameWorld'',' num2str(ndx) ');'],'interpreter','none')
                set(xl,'buttondownfcn',['virmenEventHandler(''editWorld'',' num2str(ndx) ');'])
                
                xs(end+1) = range(get(ax,'xlim')); %#ok<AGROW>
                ys(end+1) = range(get(ax,'ylim')); %#ok<AGROW>
                zs(end+1) = range(get(ax,'zlim')); %#ok<AGROW>
                axs(end+1) = ax; %#ok<AGROW>
            else
                axis(ax,'off');
            end
            ndx = ndx+1;
        end
    end
    
    if ~isempty(xs)
        xs = max(xs)/2*1.05;
        ys = max(ys)/2*1.05;
        zs = max(zs)/2*1.05;
        for ndx = 1:length(axs);
            xl = get(axs(ndx),'xlim');
            yl = get(axs(ndx),'ylim');
            zl = get(axs(ndx),'zlim');
            set(axs(ndx),'xlim',mean(xl)+[-xs xs]);
            set(axs(ndx),'ylim',mean(yl)+[-ys ys]);
            set(axs(ndx),'zlim',mean(zl)+[-zs zs]);
        end
    end
end

f = handles.buttons.(makeVar('Show/hide triangulation'));
if handles.state.showTriangulation == 1
    set(f,'state','on');
else
    set(f,'state','off');
end

f = handles.buttons.(makeVar('Change world background'));
c = repmat(permute(handles.exper.worlds{wNum}.backgroundColor,[1 3 2]),[16 16 1]);
c([1 end],:,:) = 0;
c(:,[1 end],:) = 0;
set(f,'cdata',c);

f = handles.buttons.(makeVar('Wireframe on/off'));
if handles.state.showWireframe == 1
    set(f,'state','on');
else
    set(f,'state','off');
end

f = handles.buttons.(makeVar('World transparency on/off'));
onCall = get(f,'oncallback');
offCall = get(f,'offcallback');
set(f,'oncallback','%','offcallback','%');
if handles.exper.worlds{wNum}.transparency == 1
    set(f,'state','on');
else
    set(f,'state','off');
end
set(f,'oncallback',onCall,'offcallback',offCall);

drawnow;
if ~strcmp(handles.highlightFig.(type),'0')
    set(findobj(guifig,'type','uipanel'),'shadowcolor',[.5 .5 .5],'highlightcolor','w');
    set(handles.figs.(handles.highlightFig.(type)),'shadowcolor','r','highlightcolor','r');
end
if isempty(findobj(guifig,'type','uipanel','visible','on','shadowcolor','r'))
    priority = {'textureSketch','worldSketch','worldsMenu','textureDrawing', ...
        'worldDrawing','shapeProperties','objectProperties','experimentProperties','variablesTable'};
    set(findobj(guifig,'type','uipanel'),'shadowcolor',[.5 .5 .5],'highlightcolor','w');
    for ndx = 1:length(priority)
        if strcmp(get(handles.figs.(priority{ndx}),'visible'),'on')
            set(handles.figs.(priority{ndx}),'shadowcolor','r','highlightcolor','r');
            break
        end
    end
end

set(handles.separated,'separator','on');
bNames = fieldnames(handles.buttons);
for ndx = 1:length(bNames)
    ud = get(handles.buttons.(bNames{ndx}),'userdata');
    if ~strcmp(ud,'n/a') && ~isempty(ud)
        ud(strfind(ud,'"')) = [];
        st = strfind(ud,',');
        st = [0 st length(ud)+1]; %#ok<AGROW>
        isVisible = false;
        for fg = 1:length(st)-1
            if all(get(handles.figs.(ud(st(fg)+1:st(fg+1)-1)),'highlightcolor')==[1 0 0])
                isVisible = true;
            end
        end
        if isVisible
            set(handles.buttons.(bNames{ndx}),'visible','on');
        else
            set(handles.buttons.(bNames{ndx}),'visible','off','separator','off');
        end
    end
end

mNames = fieldnames(handles.menus);
for ndx = 1:length(mNames)
    if ~strcmp(get(handles.menus.(mNames{ndx}),'userdata'),'n/a') && ~isempty(get(handles.menus.(mNames{ndx}),'userdata'))
        f = handles.buttons.(mNames{ndx});
        if strcmp(get(f,'enable'),'on') && strcmp(get(f,'visible'),'on')
            set(handles.menus.(mNames{ndx}),'enable','on');
            if strcmp(get(f,'type'),'uitoggletool')
                set(handles.menus.(mNames{ndx}),'checked',get(f,'state'));
            end
        else
            set(handles.menus.(mNames{ndx}),'enable','off');
        end
    end
end

% Resize tables
wdt.table_variables = [.5 .5];
wdt.table_objectProperties = [.5 .5];
wdt.table_objectLocations = [.16 .42 .42];
wdt.table_shapeProperties = [.5 .5];
wdt.table_shapeLocations = [.2 .4 .4];
tabs = findobj(guifig,'type','uitable');
set(tabs,'units','pixels');
for ndx = 1:length(tabs)
    ps = get(tabs(ndx),'position');
    w = round((ps(3)-20)*wdt.(get(tabs(ndx),'tag')));
    w(w>150) = 150;
    set(tabs(ndx),'columnwidth',num2cell(w))
end
set(tabs,'units','normalized');

guidata(guifig, handles);

function str = makeVar(str)

str = str(regexp(str,'[A-Za-z]'));