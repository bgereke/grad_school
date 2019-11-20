function [figs, buttons, menus] = createFigures

% Read GUI data file
fid = fopen('guidata.txt');
txt = textscan(fid,'%s','delimiter','\t');
txt = txt{1};
numCol = find(cellfun(@(x)strcmp(x,'Icon'),txt),1);
txt = reshape(txt,numCol,length(txt)/numCol)';
fclose(fid);

% Read headers
colm = [];
for ndx = 1:size(txt,2)
    colm.(txt{1,ndx}) = ndx;
end
txt = txt(2:end,:);

global guifig;
handles = guidata(guifig);
set(guifig,'keypressfcn',@keyboardShortcuts);
toolbar = uitoolbar(guifig);

panels = findall(guifig,'type','uipanel');
for ndx = 1:length(panels)
    tag = get(panels(ndx),'tag');
    figs.(tag) = panels(ndx);
end

str = txt(:,colm.Menu);
menus = struct;
for ndx = 1:length(str)
    if ~isfield(menus,str{ndx})
        menus.(str{ndx}) = uimenu('Label',str{ndx});
    end
end

mfile = mfilename('fullpath');
path = fileparts(mfile);
mt = dir([path filesep 'icons' filesep '*.png']);
if ~exist([path filesep 'icons' filesep 'allicons.mat'],'file')
    needNewIcons = true;
    for ndx = 1:length(mt)
        allIcons.(mt(ndx).name(1:end-4)) = imread([path filesep 'icons' filesep mt(ndx).name]);
    end
else
    needNewIcons = false;
    load([path filesep 'icons' filesep 'allicons.mat']);
end

for ndx = 1:size(txt,1)
    row = txt(ndx,:);
    switch row{colm.Type}
        case 'button'
            if ~strcmp(row{colm.Icon},'n/a')
                h = uipushtool(toolbar,'tooltipstring',row{colm.ToolTip});
                if needNewIcons
                    f = strfind(row{colm.Icon},'_');
                    if ~isempty(f)
                        main = allIcons.(row{colm.Icon}(f+1:end));
                        inset = allIcons.(row{colm.Icon}(1:f-1));
                        allIcons.(row{colm.Icon}) = createCombinedIcon(main,inset);
                    end
                end
                set(h,'cdata',allIcons.(row{colm.Icon}));
                set(h,'clickedcallback',['virmenEventHandler(''' row{colm.Callback} ''',''' row{colm.CallbackArgument} ''');']);
                set(h,'separator',row{colm.Separator});
                set(h,'userdata',row{colm.Figure});
                if strcmp(row{colm.ToolTip},'Busy')
                    set(h,'visible','off');
                end
            end
            
            m = uimenu(menus.(row{colm.Menu}),'Label',row{colm.MenuLabel},'callback', ...
                ['virmenEventHandler(''' row{colm.Callback} ''',''' row{colm.CallbackArgument} ''');'],'userdata',row{colm.ToolTip});
            if length(get(menus.(row{colm.Menu}),'children')) > 1
                set(m,'separator',row{colm.Separator});
            end
            
            buttons.(makeVar(row{colm.ToolTip})) = h;
            menus.(makeVar(row{colm.ToolTip})) = m;
        case 'dropdown'
            hm = uisplittool(toolbar,'tooltipstring',row{colm.ToolTip});
            set(hm,'separator',row{colm.Separator});
            f = strfind(row{colm.Icon},'_');
            if isempty(f)
                set(hm,'cdata',allIcons.(row{colm.Icon}));
            else
                main = allIcons.(row{colm.Icon}(f+1:end));
                inset = allIcons.(row{colm.Icon}(1:f-1));
                set(hm,'cdata',createCombinedIcon(main,inset));
            end
            set(hm,'userdata',row{colm.Figure});
            lst = populateList(row{colm.Callback});
            drawnow
            j = get(hm,'javacontainer');
            jmenu = get(j,'menucomponent');
            h = {};
            warning off MATLAB:hg:JavaSetHGProperty;
            warning off MATLAB:hg:PossibleDeprecatedJavaSetHGProperty;
            for m = 1:length(lst)
                h{m} = jmenu.add(lst(m).name); %#ok<AGROW>
                jButton = handle(h{m}, 'CallbackProperties');
                set(jButton,'ActionPerformedCallback',['virmenEventHandler(''' lst(m).callback ''',''' lst(m).callbackArgument ''');']);
            end
            drawnow
            
            mfile = mfilename('fullpath');
            path = fileparts(mfile);
            for m = 1:length(lst)
                if ischar(lst(m).icon)
                    icon = javax.swing.ImageIcon([path filesep 'icons' filesep lst(m).icon '.png']);
                else
                    tmp = tempname;
                    imwrite(lst(m).icon,[tmp '.png']);
                    icon = javax.swing.ImageIcon([tmp '.png']);
                end
                h{m}.setIcon(icon);
            end
            
            m = uimenu(menus.(row{colm.Menu}),'Label',row{colm.MenuLabel},'userdata',row{colm.ToolTip});
            if length(get(menus.(row{colm.Menu}),'children')) > 1
                set(m,'separator',get(hm,'separator'));
            end
            for mndx = 1:length(h)
                uimenu(m,'Label',lst(mndx).name,'callback',['virmenEventHandler(''' lst(mndx).callback ''',''' lst(mndx).callbackArgument ''');'], ...
                    'accelerator',lst(mndx).shortcut,'userdata','n/a');
            end
            
            buttons.(makeVar(row{colm.ToolTip})) = hm;
            menus.(makeVar(row{colm.ToolTip})) = m;
        case 'toggle'
            if ~strcmp(row{colm.Icon},'n/a')
                h = uitoggletool(toolbar,'tooltipstring',row{colm.ToolTip});
                set(h,'separator',row{colm.Separator});
                set(h,'cdata',allIcons.(row{colm.Icon}));
                set(h,'oncallback',['virmenEventHandler(''' row{colm.Callback} ''',{''on'',''' row{colm.CallbackArgument} '''});']);
                set(h,'offcallback',['virmenEventHandler(''' row{colm.Callback} ''',{''off'',''' row{colm.CallbackArgument} '''});']);
                set(h,'userdata',row{colm.Figure});
            end
            
            m = uimenu(menus.(row{colm.Menu}),'Label',row{colm.MenuLabel},'userdata',row{colm.ToolTip}, ...
                'callback',['virmenEventHandler(''' row{colm.Callback} ''',{''switch'',''' row{colm.CallbackArgument} '''});']);
            if length(get(menus.(row{colm.Menu}),'children')) > 1
                set(m,'separator',row{colm.Separator});
            end
            
            buttons.(makeVar(row{colm.ToolTip})) = h;
            menus.(makeVar(row{colm.ToolTip})) = m;
    end
end

if needNewIcons
    save([path filesep 'icons' filesep 'allicons.mat'],'allIcons');
end

f = fieldnames(menus);
for ndx = length(f):-1:1
    if strcmp(get(menus.(f{ndx}),'label'),'NA')
        delete(menus.(f{ndx}));
    end
end
for s = 1:length(handles.shortcuts)
    if strcmp(handles.shortcuts(s).modifier,'control') && length(handles.shortcuts(s).key)==1
        for ndx = 1:length(f)
            if strcmp(get(menus.(f{ndx}),'callback'),['virmenEventHandler(''' handles.shortcuts(s).callback ''',''' handles.shortcuts(s).input ''');']) || ...
                    strcmp(get(menus.(f{ndx}),'callback'),['virmenEventHandler(''' handles.shortcuts(s).callback ''',{''switch'',''' handles.shortcuts(s).input '''});'])
                set(menus.(f{ndx}),'accelerator',handles.shortcuts(s).key);
            end
        end
    end
end

% Layout all the figures
layout = layoutList;
layout = layout{1};
layout = rmfield(layout,'name');
layout = rmfield(layout,'icon');
figureLayout(figs,layout);

f = menus.(makeVar('Layout'));
ch = get(f,'children');
set(ch(end),'checked','on');
f = menus.(makeVar('Experiment layout'));
set(f,'checked','on');

guidata(guifig, handles)

function icon = createCombinedIcon(main,inset)

icon = uint8(240*ones(16,16,3));
main = imresize(main,[12 12]);
icon(1:12,end-11:end,:) = main;
inset = imresize(inset,[10 10]);
icon(end-9:end,1:10,:) = inset;

function str = makeVar(str)

str = str(regexp(str,'[A-Za-z]'));