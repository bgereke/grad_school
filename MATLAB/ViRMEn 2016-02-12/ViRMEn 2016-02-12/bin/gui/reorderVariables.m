function varargout = reorderVariables(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @reorderVariables_OpeningFcn, ...
                   'gui_OutputFcn',  @reorderVariables_OutputFcn, ...
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


% --- Executes just before reorderVariables is made visible.
function reorderVariables_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>

backupUnits = get(0,'units');
set(0,'units','pixels');
scr = get(0,'screensize');
pos = get(handles.mainFigure,'position');
set(handles.mainFigure,'position',[(scr(3:4)-pos(3:4))/2 pos(3) pos(4)]);
set(0,'units',backupUnits);

set(handles.list_variables,'string',varargin{1},'value',1);
set(handles.push_up,'enable','off');
handles.order = 1:length(varargin{1});

set(handles.text_Label,'string',varargin{2});
set(handles.mainFigure,'name',['Sort ' lower(varargin{2})]);

mfile = mfilename('fullpath');
mfile = fileparts(mfile);
if ismac
    mfile = [filesep mfile];
end
path = [mfile filesep 'icons' filesep 'moveUp.png'];
iconUrl = strrep(['file:/' path],filesep,'/');
set(handles.push_up,'string',['<html><img src="' iconUrl '"/></html>']);
path = [mfile filesep 'icons' filesep 'moveDown.png'];
iconUrl = strrep(['file:/' path],filesep,'/');
set(handles.push_down,'string',['<html><img src="' iconUrl '"/></html>']);
path = [mfile filesep 'icons' filesep 'sortAscend.png'];
iconUrl = strrep(['file:/' path],filesep,'/');
set(handles.push_AZ,'string',['<html><img src="' iconUrl '"/></html>']);
path = [mfile filesep 'icons' filesep 'sortDescend.png'];
iconUrl = strrep(['file:/' path],filesep,'/');
set(handles.push_ZA,'string',['<html><img src="' iconUrl '"/></html>']);

% Update handles structure
guidata(hObject, handles);

uiwait(handles.mainFigure);


% --- Outputs from this function are returned to the command line.
function varargout = reorderVariables_OutputFcn(hObject, eventdata, handles)  %#ok<*INUSD>

varargout{1} = handles.output;
delete(handles.mainFigure);


% --- Executes on selection change in list_variables.
function list_variables_Callback(hObject, eventdata, handles) %#ok<*DEFNU>

str = get(hObject,'string');
val = get(hObject','value');
if val == 1
    set(handles.push_up,'enable','off');
else
    set(handles.push_up,'enable','on');
end
if val == length(str)
    set(handles.push_down,'enable','off');
else
    set(handles.push_down,'enable','on');
end


% --- Executes during object creation, after setting all properties.
function list_variables_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_up.
function push_up_Callback(hObject, eventdata, handles)

str = get(handles.list_variables,'string');
val = get(handles.list_variables,'value');
str = str([1:val-2 val val-1 val+1:end]);
handles.order = handles.order([1:val-2 val val-1 val+1:end]);
set(handles.list_variables,'string',str,'value',val-1);
list_variables_Callback(handles.list_variables, eventdata, handles)
guidata(hObject, handles);

% --- Executes on button press in push_down.
function push_down_Callback(hObject, eventdata, handles)

str = get(handles.list_variables,'string');
val = get(handles.list_variables,'value');
str = str([1:val-1 val+1 val val+2:end]);
handles.order = handles.order([1:val-1 val+1 val val+2:end]);
set(handles.list_variables,'string',str,'value',val+1);
list_variables_Callback(handles.list_variables, eventdata, handles)
guidata(hObject, handles);


% --- Executes on button press in push_ok.
function push_ok_Callback(hObject, eventdata, handles)

handles.output = handles.order;
guidata(gcf, handles);
set(gcf,'closerequestfcn','%');
uiresume;

function closeFigure(hObject, eventdata, handles)

str = get(handles.list_variables,'string');
handles.output = 1:length(str);
guidata(gcf, handles);
set(gcf,'closerequestfcn','%');
uiresume;


% --- Executes on button press in push_AZ.
function push_AZ_Callback(hObject, eventdata, handles)

str = get(handles.list_variables,'string');
val = get(handles.list_variables,'value');
[str, ord] = sort(str);
handles.order = handles.order(ord);
set(handles.list_variables,'string',str,'value',find(ord==val,1));
list_variables_Callback(handles.list_variables, eventdata, handles)
guidata(hObject, handles);


% --- Executes on button press in push_ZA.
function push_ZA_Callback(hObject, eventdata, handles)

str = get(handles.list_variables,'string');
val = get(handles.list_variables,'value');
[str, ord] = sort(str);
str = str(end:-1:1);
ord = ord(end:-1:1);
handles.order = handles.order(ord);
set(handles.list_variables,'string',str,'value',find(ord==val,1));
list_variables_Callback(handles.list_variables, eventdata, handles)
guidata(hObject, handles);


% --- Executes on button press in push_cancel.
function push_cancel_Callback(hObject, eventdata, handles)

str = get(handles.list_variables,'string');
handles.output = 1:length(str);
guidata(gcf, handles);
set(gcf,'closerequestfcn','%');
uiresume;
