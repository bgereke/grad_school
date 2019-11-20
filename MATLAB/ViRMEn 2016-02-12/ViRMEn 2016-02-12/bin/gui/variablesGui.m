function varargout = variablesGui(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @variablesGui_OpeningFcn, ...
                   'gui_OutputFcn',  @variablesGui_OutputFcn, ...
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


% --- Executes just before variablesGui is made visible.
function variablesGui_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = variablesGui_OutputFcn(hObject, eventdata, handles) %#ok<*INUSD>

varargout = {};

function variablesEdit_Callback(hObject, eventdata, handles) %#ok<*DEFNU>

data = get(handles.table_variables,'data');
row = data(:,1);
indx = eventdata.Indices(1);
val  = eventdata.NewData;
global guifig;
fig = gcf;
set(0,'currentfigure',guifig);
virmenEventHandler('changeVariables',{row{indx},val});
set(0,'currentfigure',fig);
