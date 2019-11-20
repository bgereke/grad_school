function highlightFigure

global guifig
handles = guidata(guifig);

scr = get(0,'screensize');
ptr = get(0,'pointerlocation');
figpos = get(guifig,'position').*scr([3 4 3 4]);
ptr = (ptr-figpos(1:2))./figpos(3:4);
isin = structfun(@(x)all(sum(sign(cumsum(reshape(get(x,'position'),2,2),2)-[ptr' ptr']),2)'==0),handles.figs);
isvis = structfun(@(x)strcmp(get(x,'visible'),'on'),handles.figs);
f = find(and(isin,isvis),1);
if isempty(f)
    return
end
fld = fieldnames(handles.figs);
set(findobj(guifig,'type','uipanel'),'shadowcolor',[.5 .5 .5],'highlightcolor','w');
set(handles.figs.(fld{f}),'shadowcolor','r','highlightcolor','r');

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