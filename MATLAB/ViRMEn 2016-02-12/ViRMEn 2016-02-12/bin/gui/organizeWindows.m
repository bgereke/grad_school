function varargout = organizeWindows(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @organizeWindows_OpeningFcn, ...
                   'gui_OutputFcn',  @organizeWindows_OutputFcn, ...
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


% --- Executes just before organizeWindows is made visible.
function organizeWindows_OpeningFcn(hObject, ~, handles, varargin)

backupUnits = get(0,'units');
set(0,'units','pixels');
scr = get(0,'screensize');
pos = get(handles.figure_Main,'position');
set(handles.figure_Main,'position',[(scr(3:4)-pos(3:4))/2 pos(3) pos(4)]);
set(0,'units',backupUnits);

handles.windows = varargin{1};
handles.defaultProperties = varargin{2};

mfile = mfilename('fullpath');
mfile = fileparts(mfile);
if ismac
    mfile = [filesep mfile];
end
path = [mfile filesep 'icons' filesep 'monitorAdd.png'];
iconUrl = strrep(['file:/' path],filesep,'/');
set(handles.push_add,'string',['<html><img src="' iconUrl '"/></html>']);
path = [mfile filesep 'icons' filesep 'monitorDelete.png'];
iconUrl = strrep(['file:/' path],filesep,'/');
set(handles.push_delete,'string',['<html><img src="' iconUrl '"/></html>']);
path = [mfile filesep 'icons' filesep 'moveUp.png'];
iconUrl = strrep(['file:/' path],filesep,'/');
set(handles.push_up,'string',['<html><img src="' iconUrl '"/></html>']);
path = [mfile filesep 'icons' filesep 'moveDown.png'];
iconUrl = strrep(['file:/' path],filesep,'/');
set(handles.push_down,'string',['<html><img src="' iconUrl '"/></html>']);

columns = {'3D graphics','Function #','Main monitor','Monitor #','Full screen','Left','Bottom','Width','Height','Antialiasing'};
columnformat = {'logical','short','logical','short','logical','short','short','short','short','short'};
columnwidth = {'auto', 'auto', 'auto', 'auto', 'auto', 50, 50, 50, 50, 'auto'};
columneditable = true(1,10);
set(handles.table_Main,'columnname',columns,'columnformat',columnformat,'columnwidth',columnwidth,'columneditable',columneditable);

handles = updateTable(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes organizeWindows wait for user response (see UIRESUME)
uiwait(handles.figure_Main);

function handles = updateTable(handles)

data = cell(length(handles.windows),10);
for ndx = 1:length(handles.windows)
    data{ndx,1} = handles.windows{ndx}.rendering3D;
    data{ndx,2} = forceNaN(handles.windows{ndx}.transformation,~handles.windows{ndx}.rendering3D);
    data{ndx,3} = handles.windows{ndx}.primaryMonitor;
    data{ndx,4} = forceNaN(handles.windows{ndx}.monitor,handles.windows{ndx}.primaryMonitor);
    data{ndx,5} = handles.windows{ndx}.fullScreen;
    data{ndx,6} = forceNaN(handles.windows{ndx}.left,handles.windows{ndx}.fullScreen);
    data{ndx,7} = forceNaN(handles.windows{ndx}.bottom,handles.windows{ndx}.fullScreen);
    data{ndx,8} = forceNaN(handles.windows{ndx}.width,handles.windows{ndx}.fullScreen);
    data{ndx,9} = forceNaN(handles.windows{ndx}.height,handles.windows{ndx}.fullScreen);
    data{ndx,10} = handles.windows{ndx}.antialiasing;
end
set(handles.table_Main,'data',data);

function val = forceNaN(val,forcer)

if forcer
    val = NaN;
end


function editTable(hObject, eventdata, handles)  %#ok<DEFNU>

indx = eventdata.Indices(1);
switch eventdata.Indices(2)
    case 1
        handles.windows{indx}.rendering3D = eventdata.NewData;
    case 2
        handles.windows{indx}.rendering3D = true;
        handles.windows{indx}.transformation = eventdata.NewData;
    case 3
        handles.windows{indx}.primaryMonitor = eventdata.NewData;
    case 4
        handles.windows{indx}.primaryMonitor = false;
        handles.windows{indx}.monitor = eventdata.NewData;
    case 5
        handles.windows{indx}.fullScreen = eventdata.NewData;
    case 6
        handles.windows{indx}.fullScreen = false;
        handles.windows{indx}.left = eventdata.NewData;
    case 7
        handles.windows{indx}.fullScreen = false;
        handles.windows{indx}.bottom = eventdata.NewData;
    case 8
        handles.windows{indx}.fullScreen = false;
        handles.windows{indx}.width = eventdata.NewData;
    case 9
        handles.windows{indx}.fullScreen = false;
        handles.windows{indx}.height = eventdata.NewData;
    case 10
        handles.windows{indx}.antialiasing = eventdata.NewData;
end

handles = updateTable(handles);
set(handles.table_Main,'userdata',eventdata.Indices);

% Update handles structure
guidata(hObject, handles);

function selectTable(hObject, eventdata, handles)  %#ok<DEFNU>

set(hObject,'userdata',eventdata.Indices);

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = organizeWindows_OutputFcn(~, ~, handles) 

% Get default command line output from handles structure
if isempty(handles)
    varargout{1} = [];
else
    varargout{1} = handles.output;
    close(handles.figure_Main);
end


% --- Executes on button press in push_OK.
function push_OK_Callback(hObject, ~, ~) %#ok<DEFNU>

handles = guidata(hObject);
handles.output = handles.windows;
guidata(hObject, handles);
uiresume;


% --- Executes on button press in push_add.
function push_add_Callback(hObject, ~, handles) %#ok<DEFNU>

handles.windows{end+1} = virmenWindow;
indx = length(handles.windows);
handles.windows{end}.rendering3D = handles.defaultProperties.windowRendering3D(min(end,indx));
handles.windows{end}.transformation = handles.defaultProperties.windowTransformation(min(end,indx));
handles.windows{end}.primaryMonitor = handles.defaultProperties.windowPrimaryMonitor(min(end,indx));
handles.windows{end}.monitor = handles.defaultProperties.windowMonitor(min(end,indx));
handles.windows{end}.fullScreen = handles.defaultProperties.windowFullScreen(min(end,indx));
handles.windows{end}.left = handles.defaultProperties.windowLeft(min(end,indx));
handles.windows{end}.bottom = handles.defaultProperties.windowBottom(min(end,indx));
handles.windows{end}.width = handles.defaultProperties.windowWidth(min(end,indx));
handles.windows{end}.height = handles.defaultProperties.windowHeight(min(end,indx));
handles.windows{end}.antialiasing = handles.defaultProperties.windowAntialiasing(min(end,indx));

handles = updateTable(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in push_delete.
function push_delete_Callback(hObject, ~, handles) %#ok<DEFNU>

pos = get(handles.table_Main,'userdata');
if isempty(pos)
    errordlg('No window selected.','Error');
    return
end
if length(handles.windows) == 1
    errordlg('At least one window is required.','Error');
    return
end
handles.windows(pos(1)) = [];

handles = updateTable(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in push_cancel.
function push_cancel_Callback(hObject, ~, ~) %#ok<DEFNU>

handles = guidata(hObject);
handles.output = [];
guidata(hObject, handles);
uiresume;


% --- Executes on button press in push_up.
function push_up_Callback(hObject, ~, handles) %#ok<DEFNU>

pos = get(handles.table_Main,'userdata');
if isempty(pos)
    return
end
if pos(1) == 1
    return
end
handles.windows = handles.windows([1:pos(1)-2 pos(1) pos(1)-1 pos(1)+1:end]);
handles = updateTable(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in push_down.
function push_down_Callback(hObject, ~, handles) %#ok<DEFNU>

pos = get(handles.table_Main,'userdata');
if isempty(pos)
    return
end
if pos(1) == length(handles.windows)
    return
end
handles.windows = handles.windows([1:pos(1)-1 pos(1)+1 pos(1) pos(1)+2:end]);
handles = updateTable(handles);

% Update handles structure
guidata(hObject, handles);
