function varargout = refiningGui(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @refiningGui_OpeningFcn, ...
                   'gui_OutputFcn',  @refiningGui_OutputFcn, ...
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


% --- Executes just before refiningGui is made visible.
function refiningGui_OpeningFcn(hObject, ~, handles, varargin)

set(handles.figure_main,'name',['Refining ' varargin{1}.fullName]);
set(handles.figure_main,'closerequestfcn','handles.output = []; guidata(gcf, handles); uiresume;');

handles.texture = copyVirmenObject(varargin{1});
if ~isempty(varargin{1}.parent) && ~isempty(varargin{1}.parent.parent)
    handles.backgroundColor = varargin{1}.parent.parent.backgroundColor;
else
    handles.backgroundColor = [0 0 0];
end
handles.triangulationColor = varargin{2};
handles.showTriangulation = true;

handles.minGrid = 1;
handles.maxGrid = 100;
handles.minRefining = 1;
handles.maxRefining = 100;

set(handles.slider_verticalGrid,'min',handles.minGrid,'max',handles.maxGrid,'sliderstep',[1 5]/(handles.maxGrid-handles.minGrid));
set(handles.slider_horizontalGrid,'min',handles.minGrid,'max',handles.maxGrid,'sliderstep',[1 5]/(handles.maxGrid-handles.minGrid));
set(handles.slider_verticalRefining,'min',handles.minRefining,'max',handles.maxRefining,'sliderstep',[1 5]/(handles.maxRefining-handles.minRefining));
set(handles.slider_horizontalRefining,'min',handles.minRefining,'max',handles.maxRefining,'sliderstep',[1 5]/(handles.maxRefining-handles.minRefining));

handles = updateControls(handles);
guidata(hObject, handles);

% UIWAIT makes refiningGui wait for user response (see UIRESUME)
uiwait(handles.figure_main);

function handles = updateControls(handles)

set(handles.check_verticalTilable,'value',handles.texture.tilable(1));
set(handles.check_horizontalTilable,'value',handles.texture.tilable(2));
set(handles.edit_verticalRefining,'string',handles.texture.getValue.refining{1});
set(handles.edit_horizontalRefining,'string',handles.texture.getValue.refining{2});
set(handles.edit_verticalGrid,'string',handles.texture.getValue.grid{1});
set(handles.edit_horizontalGrid,'string',handles.texture.getValue.grid{2});
set(handles.slider_verticalRefining,'value',max(handles.minRefining,min(handles.maxRefining,1/handles.texture.refining(1))));
set(handles.slider_horizontalRefining,'value',max(handles.minRefining,min(handles.maxRefining,1/handles.texture.refining(2))));
set(handles.slider_verticalGrid,'value',max(handles.minGrid,min(handles.maxGrid,handles.texture.grid(1))));
set(handles.slider_horizontalGrid,'value',max(handles.minGrid,min(handles.maxGrid,handles.texture.grid(2))));

handles.texture.compute;

subplot(handles.axes_texture)
h = handles.texture.draw;
set(h,'edgecolor','r');
view(2);
xl = [-0.025*handles.texture.width 1.025*handles.texture.width];
yl = [-0.025*handles.texture.height 1.025*handles.texture.height];
set(gca,'units','pixels');
pos = get(gca,'position');
set(gca,'units','normalized');
ar = pos(4)/pos(3);
if range(yl)/range(xl) > ar
    xl = mean(xl) + (xl-mean(xl))/range(xl)*range(yl)/ar;
else
    yl = mean(yl) + (yl-mean(yl))/range(yl)*range(xl)*ar;
end
xlim(xl);
ylim(yl);
set(gca,'color',handles.backgroundColor);
set(gca,'xtick',[],'ytick',[],'box','on');

num = length(find(~isnan(handles.texture.triangles.cdata(handles.texture.triangles.triangulation(:,1),1))));
title([num2str(num) ' visible triangles']);
if handles.showTriangulation == 0
    set(h,'edgecolor','none');
else
    set(h,'edgecolor',handles.triangulationColor);
end


% --- Outputs from this function are returned to the command line.
function varargout = refiningGui_OutputFcn(~, ~, handles) 

varargout{1} = handles.output;
set(gcf,'closerequestfcn','%');
delete(gcf);


function slider_verticalRefining_Callback(hObject, ~, handles)

handles.texture.refining(1) = 1/get(hObject,'value');
handles = updateControls(handles);
guidata(hObject, handles);


function slider_horizontalRefining_Callback(hObject, ~, handles)

handles.texture.refining(2) = 1/get(hObject,'value');
handles = updateControls(handles);
guidata(hObject, handles);

function slider_verticalGrid_Callback(hObject, ~, handles)

handles.texture.grid(1) = round(get(hObject,'value'));
handles = updateControls(handles);
guidata(hObject, handles);

function slider_horizontalGrid_Callback(hObject, ~, handles)

handles.texture.grid(2) = round(get(hObject,'value'));
handles = updateControls(handles);
guidata(hObject, handles);


function edit_verticalRefining_Callback(hObject, ~, handles)

val = handles.texture.getValue.refining;
val{1} = get(hObject,'string');
handles.texture.refining = val;
handles = updateControls(handles);
guidata(hObject, handles);

function edit_horizontalRefining_Callback(hObject, ~, handles)

val = handles.texture.getValue.refining;
val{2} = get(hObject,'string');
handles.texture.refining = val;
handles = updateControls(handles);
guidata(hObject, handles);

function edit_verticalGrid_Callback(hObject, ~, handles)

val = handles.texture.getValue.grid;
val{1} = get(hObject,'string');
handles.texture.grid = val;
handles = updateControls(handles);
guidata(hObject, handles);

function edit_horizontalGrid_Callback(hObject, ~, handles)

val = handles.texture.getValue.grid;
val{2} = get(hObject,'string');
handles.texture.grid = val;
handles = updateControls(handles);
guidata(hObject, handles);

function check_verticalTilable_Callback(hObject, ~, handles)

handles.texture.tilable(1) = get(hObject,'value');
handles = updateControls(handles);
guidata(hObject, handles);

function check_horizontalTilable_Callback(hObject, ~, handles)

handles.texture.tilable(2) = get(hObject,'value');
handles = updateControls(handles);
guidata(hObject, handles);

function push_ok_Callback(hObject, ~, handles)

handles.output = {get(handles.edit_verticalRefining,'string') ...
    get(handles.edit_horizontalRefining,'string') ...
    get(handles.edit_verticalGrid,'string') ...
    get(handles.edit_horizontalGrid,'string') ...
    num2str(get(handles.check_verticalTilable,'value')) ...
    num2str(get(handles.check_horizontalTilable,'value'))};
guidata(hObject, handles);
uiresume;

function push_cancel_Callback(hObject, ~, handles)

handles.output = [];
guidata(hObject, handles);
uiresume;







% --- Executes during object creation, after setting all properties.
function slider_verticalRefining_CreateFcn(hObject, ~, ~) %#ok<*DEFNU>
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function slider_verticalGrid_CreateFcn(hObject, ~, ~)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function edit_verticalRefining_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_verticalGrid_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_horizontalGrid_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_horizontalRefining_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function slider_horizontalGrid_CreateFcn(hObject, ~, ~)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function slider_horizontalRefining_CreateFcn(hObject, ~, ~)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
