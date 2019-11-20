function varargout = chooseObjectGui(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @chooseObjectGui_OpeningFcn, ...
                   'gui_OutputFcn',  @chooseObjectGui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before chooseObjectGui is made visible.
function chooseObjectGui_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>

handles.currentExperiment = varargin{1};
handles.type = varargin{2};

backupUnits = get(0,'units');
set(0,'units','pixels');
scr = get(0,'screensize');
set(handles.mainFigure,'position',[50 100 scr(3)-100 scr(4)-200]);
set(0,'units',backupUnits);

mfile = mfilename('fullpath');
path = fileparts(mfile);
mt = dir([path filesep '..' filesep '..' filesep 'experiments' filesep '*.mat']);

str = {'[Current experiment]'};
for ndx = 1:length(mt)
    str{ndx+1} = mt(ndx).name(1:end-4);
end
set(handles.list_experiments,'string',str);

handles = loadExperiment(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes chooseObjectGui wait for user response (see UIRESUME)
uiwait(handles.mainFigure);


function handles = loadExperiment(handles)

val = get(handles.list_experiments,'value');
if val == 1
    exper = handles.currentExperiment;
else
    str = get(handles.list_experiments,'string');
    load(str{val});
end

items = exper.descendants;
for ndx = length(items):-1:1
    supr = superclasses(items{ndx});
    if length(setdiff(supr,handles.type))==length(supr) && ~strcmp(class(items{ndx}),handles.type)
        items(ndx) = [];
    end
end

switch handles.type
    case 'virmenTexture'
        for ndx2 = length(items):-1:1
            for ndx1 = 1:ndx2-1
                if isequalwithequalnans(items{ndx2}.triangles,items{ndx1}.triangles) %#ok<DISEQN>
                    items(ndx2) = [];
                    break
                end
            end
        end
end

set(handles.mainFigure,'userdata',items);
handles.page = 1;
handles = drawPage(handles);

function handles = drawPage(handles)

items = get(handles.mainFigure,'userdata');

switch handles.type
    case 'virmenWorld'
        wd = 3;
        hg = 2;
    otherwise
        wd = 4;
        hg = 3;
end

if handles.page < 1
    handles.page = ceil(length(items)/(wd*hg));
end
if (handles.page-1)*wd*hg > length(items)
    handles.page = 1;
end
set(handles.text_pageNum,'string',[num2str(handles.page) ' / ' num2str(ceil(length(items)/(wd*hg)))]);

scale = .95;

listPos = get(handles.list_experiments,'position');
ndx = (wd*hg)*(handles.page-1)+1;
delete(findobj(handles.mainFigure,'type','axes'));
xs = [];
ys = [];
zs = [];
axs = [];
for y = hg:-1:1
    for x = 1:wd
        pos = [(x-1)/wd+(1-scale)/(2*wd) (y-1)/hg+(1-scale)/(2*hg) scale/wd scale/hg];
        pos(1) = (listPos(1)+listPos(3))+pos(1)*(1-(listPos(1)+listPos(3)));
        pos(2) = listPos(2)+pos(2)*listPos(4);
        pos(3) = pos(3)*(1-(listPos(1)+listPos(3)));
        pos(4) = pos(4)*listPos(4);
        axes('outerposition',pos);
        if ndx <= length(items)
            switch handles.type
                case 'virmenWorld'
                    [h, he, hp] = items{ndx}.draw2D;
                    if ~isempty(h)
                        delete([he hp]);
                        he = [];
                    end
                    view(3)
                    axis tight;
                    axis equal;
                    set(gca,'xtick',[],'ytick',[],'ztick',[],'box','off');
                    tt = title(items{ndx}.fullName);
                    hands = [h he tt];
                    
                    xs(end+1) = range(get(gca,'xlim')); %#ok<AGROW>
                    ys(end+1) = range(get(gca,'ylim')); %#ok<AGROW>
                    zs(end+1) = range(get(gca,'zlim')); %#ok<AGROW>
                    axs(end+1) = gca; %#ok<AGROW>
                case 'virmenTexture'
                    h = items{ndx}.draw;
                    set(h,'edgecolor','none');
                    axis tight;
                    xl = xlim;
                    yl = ylim;
                    axis auto;
                    mx = max([range(xl) range(yl)]);
                    view(2);
                    axis off;
                    xlim([mean(xl)-mx/2 mean(xl)+mx/2]);
                    ylim([mean(yl)-mx/2 mean(yl)+mx/2]);
                    axis equal;
                    tt = title(items{ndx}.fullName);
                    hands = [h tt];
                case 'virmenObject'
                    h = items{ndx}.draw3D;
                    set(h,'edgecolor','none');
                    axis tight;
                    axis equal;
                    set(gca,'xtick',[],'ytick',[],'ztick',[]);
                    grid off;
                    set(gca,'color','none');
                    tt = title(items{ndx}.fullName);
                    hands = [h tt];
            end
            set(tt,'interpreter','none','units','pixels');
            hands = [hands gca]; %#ok<AGROW>
            set(hands,'userdata',ndx);
            set(hands,'buttondownfcn',@clickObject);
        else
            axis off
        end
        ndx = ndx+1;
    end
end

switch handles.type
    case 'virmenWorld'
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

function clickObject(src,evt)

handles = guidata(src);
handles.output = get(src,'userdata');
guidata(src, handles);
uiresume;

% --- Outputs from this function are returned to the command line.
function varargout = chooseObjectGui_OutputFcn(hObject, eventdata, handles)  %#ok<*INUSD>

if isempty(handles)
    varargout{1} = [];
    varargout{2} = [];
else
    items = get(gcf,'userdata');
    obj = items{handles.output};
    
    vars = variablesList(obj);
    
    vName = fieldnames(vars);
    for ndx = 1:length(vName)
        if isfield(handles.currentExperiment.variables,vName{ndx})
            if ~strcmp(handles.currentExperiment.variables.(vName{ndx}),vars.(vName{ndx}))
                str1 = ['Old: ' handles.currentExperiment.variables.(vName{ndx})];
                str2 = ['New: ' vars.(vName{ndx})];
                button = questdlg(['Which value of variable ' vName{ndx} ' should be assigned?'],'Variable import',str1,str2,str1);
                if isempty(button)
                    varargout{1} = [];
                    varargout{2} = [];
                    return
                end
                if strcmp(button,str1)
                    vars = rmfield(vars,vName{ndx});
                end
            end
        end
    end
    
    varargout{1} = obj;
    varargout{2} = vars;
    close(gcf)
end


% --- Executes on selection change in list_experiments.
function list_experiments_Callback(hObject, eventdata, handles) %#ok<*DEFNU>

handles = loadExperiment(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function list_experiments_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_left.
function push_left_Callback(hObject, eventdata, handles)

handles.page = handles.page - 1;
handles = drawPage(handles);
guidata(hObject, handles);


% --- Executes on button press in push_right.
function push_right_Callback(hObject, eventdata, handles)

handles.page = handles.page + 1;
handles = drawPage(handles);
guidata(hObject, handles);
