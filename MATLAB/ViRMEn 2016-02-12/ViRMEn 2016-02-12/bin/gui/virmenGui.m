function varargout = virmenGui(varargin)
% *************************************************************************
% Copyright 2013, Princeton University.  All rights reserved.
% 
% By using this software the USER indicates that he or she has read,
% understood and will comply with the following:
% 
%  --- Princeton University hereby grants USER nonexclusive permission to
% use, copy and/or modify this software for internal, noncommercial,
% research purposes only. Any distribution, including publication or
% commercial sale or license, of this software, copies of the software, its
% associated documentation and/or modifications of either is strictly
% prohibited without the prior consent of Princeton University. Title to
% copyright to this software and its associated documentation shall at all
% times remain with Princeton University.  Appropriate copyright notice
% shall be placed on all software copies, and a complete copy of this
% notice shall be included in all copies of the associated documentation.
% No right is granted to use in advertising, publicity or otherwise any
% trademark, service mark, or the name of Princeton University. 
% 
%  --- This software and any associated documentation is provided "as is" 
% 
% PRINCETON UNIVERSITY MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR
% IMPLIED, INCLUDING THOSE OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR
% PURPOSE, OR THAT  USE OF THE SOFTWARE, MODIFICATIONS, OR ASSOCIATED
% DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS OR
% OTHER INTELLECTUAL PROPERTY RIGHTS OF A THIRD PARTY. 
% 
% Princeton University shall not be liable under any circumstances for any
% direct, indirect, special, incidental, or consequential damages with
% respect to any claim by USER or any third party on account of or arising
% from the use, or inability to use, this software or its associated
% documentation, even if Princeton University has been advised of the
% possibility of those damages.
% *************************************************************************

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @virmenGui_OpeningFcn, ...
                   'gui_OutputFcn',  @virmenGui_OutputFcn, ...
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


% --- Executes just before virmenGui is made visible.
function virmenGui_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>

global guifig

guifig = findall(0,'name','ViRMEn');

set(guifig,'position',[.1 .15 .8 .75]);

% Read shortcuts
fid = fopen('keyboardShortcuts.txt');
txt = textscan(fid,'%s','delimiter','\t');
txt = txt{1};
numCol = 5;
txt = reshape(txt,numCol,length(txt)/numCol)';
fclose(fid);
handles.shortcuts = struct;
for ndx = 2:size(txt,1)
    for j = 1:size(txt,2)
        handles.shortcuts(ndx-1).(txt{1,j}) = txt{ndx,j};
    end
end

% Read event booleans file
fid = fopen('eventBooleans.txt');
txt = textscan(fid,'%s','delimiter','\t');
txt = txt{1};
indx = find(cellfun(@(x)isnan(str2double(x)),txt));
numCol = length(txt)-indx(end)+1;
txt = reshape(txt,numCol,length(txt)/numCol)';
for ndx = 2:size(txt,1)
    handles.evtNames.(txt{ndx,1}) = ndx-1;
    handles.historyBool.(txt{ndx,1}) = str2double(txt{ndx,2});
    handles.highlightFig.(txt{ndx,1}) = txt{ndx,3};
end
for ndx = 4:size(txt,2)
    handles.figNames.(txt{1,ndx}) = ndx-3;
end
handles.bools = cellfun(@str2double,txt(2:end,4:end));
handles.changedFigs = zeros(1,size(handles.bools,2));
fclose(fid);

% Create macros structure
handles.macros = struct;

mfile = mfilename('fullpath');
path = fileparts(mfile);
if ismac && ~exist('/usr/local/lib/libglfw.3.dylib','file')
    try
        myCopyFile([path '/../engine/libglfw.3.dylib'], '/usr/local/lib/libglfw.3.dylib');
    catch %#ok<CTCH>
        msg = ['To run ViRMEn you need to copy the GLFW graphics library ' ...
               'to your local libraries directory, which requires ' ...
               'administrator priveledges. Press OK to continue, and you will be' ...
               'prompted for your administrator password. Press Cancel to exit ViRMEn.'];
        choice = questdlg(msg, 'ViRMEn Installation', 'OK', 'Cancel', 'OK');
        switch choice
            case 'OK'
                eval(['!sudo cp ' path '/../engine/libglfw.3.dylib /usr/local/lib']);
            case 'Cancel'
                disp('ViRMEn installation canceled.');
                return;
        end
    end
end

foldNames = {'defaults','experiments','macros','movements','objects','shapes','transformations'};
for fNdx = 1:length(foldNames)
    if ~exist([path filesep '..' filesep '..' filesep foldNames{fNdx}],'dir')
        mkdir([path filesep '..' filesep '..' filesep foldNames{fNdx}]);
    end
end

if ~exist([path filesep '..' filesep '..' filesep 'defaults' filesep 'defaultVirmenCode.m'],'file')
    myCopyFile([path filesep 'defaultVirmenCode.m'],[path filesep '..' filesep '..' filesep 'defaults' filesep 'defaultVirmenCode.m']);
end
if ~exist([path filesep '..' filesep '..' filesep 'defaults' filesep 'defaultLayout.txt'],'file')
    myCopyFile([path filesep 'defaultLayout.txt'],[path filesep '..' filesep '..' filesep 'defaults' filesep 'defaultLayout.txt']);
end
if ~exist([path filesep '..' filesep '..' filesep 'defaults' filesep 'defaultFunctions.txt'],'file')
    myCopyFile([path filesep 'defaultFunctions.txt'],[path filesep '..' filesep '..' filesep 'defaults' filesep 'defaultFunctions.txt']);
end
if ~exist([path filesep '..' filesep '..' filesep 'defaults' filesep 'defaultProperties.txt'],'file')
    myCopyFile([path filesep 'defaultProperties.txt'],[path filesep '..' filesep '..' filesep 'defaults' filesep 'defaultProperties.txt']);
end

% Transfer built-on objects and shapes
mt = dir([path filesep '..' filesep 'gui' filesep 'builtinShapes' filesep 'shape*.*']);
for ndx = 1:length(mt)
    toQuit = manageFile([path filesep '..' filesep 'gui' filesep 'builtinShapes' filesep mt(ndx).name], ...
        [path filesep '..' filesep '..' filesep 'shapes' filesep mt(ndx).name]);
    if toQuit
        delete(guifig);
        return;
    end
end
mt = dir([path filesep '..' filesep 'gui' filesep 'builtinObjects' filesep 'object*.*']);
for ndx = 1:length(mt)
    toQuit = manageFile([path filesep '..' filesep 'gui' filesep 'builtinObjects' filesep mt(ndx).name], ...
        [path filesep '..' filesep '..' filesep 'objects' filesep mt(ndx).name]);
    if toQuit
        delete(guifig);
        return;
    end
end
mt = dir([path filesep '..' filesep 'gui' filesep 'builtinTransformations' filesep 'transform*.*']);
for ndx = 1:length(mt)
    toQuit = manageFile([path filesep '..' filesep 'gui' filesep 'builtinTransformations' filesep mt(ndx).name], ...
        [path filesep '..' filesep '..' filesep 'transformations' filesep mt(ndx).name]);
    if toQuit
        delete(guifig);
        return;
    end
end
mt = dir([path filesep '..' filesep 'gui' filesep 'builtinMovements' filesep 'move*.*']);
for ndx = 1:length(mt)
    toQuit = manageFile([path filesep '..' filesep 'gui' filesep 'builtinMovements' filesep mt(ndx).name], ...
        [path filesep '..' filesep '..' filesep 'movements' filesep mt(ndx).name]);
    if toQuit
        delete(guifig);
        return;
    end
end

% Read default properties
fid = fopen([path filesep 'defaultProperties.txt']);
txt = textscan(fid,'%s','delimiter','\t');
txt = txt{1};
for ndx = 1:2:length(txt)-1
    handles.defaultProperties.(txt{ndx}) = eval(txt{ndx+1});
end
fid = fopen([path filesep '..' filesep '..' filesep 'defaults' filesep 'defaultProperties.txt']);
txt = textscan(fid,'%s','delimiter','\t');
txt = txt{1};
for ndx = 1:2:length(txt)-1
    handles.defaultProperties.(txt{ndx}) = eval(txt{ndx+1});
end
fclose(fid);

handles.layouts = layoutList;

handles.state = virmenGuiState(handles.defaultProperties);
handles.history.states = cell(1,100);
for ndx = 1:length(handles.history.states)
    handles.history.states{ndx}.state = [];
    handles.history.states{ndx}.undo = struct;
    handles.history.states{ndx}.redo = struct;
end
handles.history.position = 0;
handles.exper = virmenExperiment;

prop = fieldnames(handles.defaultProperties);
numWin = 1;
for ndx = 1:length(prop)
    if length(prop{ndx})>6 && strcmp(prop{ndx}(1:6),'window')
        numWin = max(numWin,length(handles.defaultProperties.(prop{ndx})));
    end
end
for ndx = 2:numWin
    handles.exper.windows{ndx} = virmenWindow;
end
for ndx = 1:length(handles.exper.windows)
    handles.exper.windows{ndx}.rendering3D = handles.defaultProperties.windowRendering3D(min(end,ndx));
    handles.exper.windows{ndx}.transformation = handles.defaultProperties.windowTransformation(min(end,ndx));
    handles.exper.windows{ndx}.primaryMonitor = handles.defaultProperties.windowPrimaryMonitor(min(end,ndx));
    handles.exper.windows{ndx}.monitor = handles.defaultProperties.windowMonitor(min(end,ndx));
    handles.exper.windows{ndx}.fullScreen = handles.defaultProperties.windowFullScreen(min(end,ndx));
    handles.exper.windows{ndx}.left = handles.defaultProperties.windowLeft(min(end,ndx));
    handles.exper.windows{ndx}.bottom = handles.defaultProperties.windowBottom(min(end,ndx));
    handles.exper.windows{ndx}.width = handles.defaultProperties.windowWidth(min(end,ndx));
    handles.exper.windows{ndx}.height = handles.defaultProperties.windowHeight(min(end,ndx));
    handles.exper.windows{ndx}.antialiasing = handles.defaultProperties.windowAntialiasing(min(end,ndx));
end
handles.exper.worlds{1}.backgroundColor = handles.defaultProperties.worldBackgroundColor;
handles.exper.worlds{1}.transparency = handles.defaultProperties.worldTransparency;
handles.exper.worlds{1}.startLocation = handles.defaultProperties.startLocation;

fid = fopen([path filesep '..' filesep '..' filesep 'defaults' filesep 'defaultFunctions.txt']);
txt = textscan(fid,'%s','delimiter','\t');
txt = txt{1};
for ndx = 1:2:length(txt)-1
    handles.exper.(txt{ndx}) = str2func(txt{ndx+1});
end
fclose(fid);

set(handles.mainFigure,'windowbuttonupfcn',@virmenDropObject,'WindowButtonMotionFcn',@virmenMoveObject);

guidata(hObject,handles);

virmenEventHandler('startProgram',[]);

function toQuit = manageFile(builtin,filename)

toQuit = false;

if ~exist(filename(1:end-8),'file')
    myCopyFile(builtin,filename(1:end-8));
    return
end
    
mtb = dir(builtin);
mtf = dir(filename(1:end-8));
if datenum(mtb.date) < datenum(mtf.date)
    [~, name, ~] = fileparts(filename);
%     errordlg(['Built in function ''' name ''' has been modified. Please rename or delete the newer user version.'],'Error','modal');
    button = questdlg(['Built-in function ''' name ''' has been modified. Choose action.'], ...
        'Question','Keep both copies','Revert to built-in copy','Quit','Keep both copies');
    switch button
        case 'Keep both copies'
            [path,file,ext] = fileparts(filename(1:end-8));
            myCopyFile(builtin,[path filesep file '_new' ext]);
            myCopyFile(builtin,filename(1:end-8));
        case 'Revert to built-in copy'
            myCopyFile(builtin,filename(1:end-8));
        case 'Quit'
            toQuit = true;
    end
else
    myCopyFile(builtin,filename(1:end-8));
    return
end


function myCopyFile(varargin)

try
    copyfile(varargin{:});
catch
    [~,filename,~] = fileparts(varargin{2});
    clear(filename);
    copyfile(varargin{:});
end

% --- Outputs from this function are returned to the command line.
function varargout = virmenGui_OutputFcn(hObject, eventdata, handles)  %#ok<*INUSD>

% Get default command line output from handles structure
varargout{1} = [];

function closereq_Callback(hObject, eventdata, handles)

virmenEventHandler('closeProgram',[]);


% --- Executes on button press in push_editCode.
function push_editCode_Callback(hObject, eventdata, handles) %#ok<*DEFNU>

virmenEventHandler('editCode','experimentCode');

% --- Executes on selection change in pop_movements.
function pop_movements_Callback(hObject, eventdata, handles)

str = get(hObject,'string');
val = get(hObject,'value');
hand = str2func(str{val});
virmenEventHandler('changeExperimentProperties',{'movementFunction',hand});


% --- Executes during object creation, after setting all properties.
function pop_movements_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_transformations.
function pop_transformations_Callback(hObject, eventdata, handles)

str = get(hObject,'string');
val = get(hObject,'value');
hand = str2func(str{val});
virmenEventHandler('changeExperimentProperties',{'transformationFunction',hand});


% --- Executes during object creation, after setting all properties.
function pop_transformations_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_object.
function pop_object_Callback(hObject, eventdata, handles)

val = get(hObject,'value');
virmenEventHandler('selectObject',{val-1,'none',1});


% --- Executes during object creation, after setting all properties.
function pop_object_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in push_objectSymbolic.
function push_objectSymbolic_Callback(hObject, eventdata, handles)


virmenEventHandler('changeSymbolicObjectLocations',[]);


function edit_tilingVertical_Callback(hObject, eventdata, handles)

virmenEventHandler('changeTiling',{1,get(hObject,'string')});


% --- Executes during object creation, after setting all properties.
function edit_tilingVertical_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_tilingHorizontal_Callback(hObject, eventdata, handles)

virmenEventHandler('changeTiling',{2,get(hObject,'string')});


% --- Executes during object creation, after setting all properties.
function edit_tilingHorizontal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_edgeRadius_Callback(hObject, eventdata, handles)

virmenEventHandler('changeObjectProperties',{'edgeRadius',get(hObject,'string')});


% --- Executes during object creation, after setting all properties.
function edit_edgeRadius_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_renameObject.
function push_renameObject_Callback(hObject, eventdata, handles)

virmenEventHandler('renameObject',[]);


function objectPropertiesEdit_Callback(hObject, eventdata, handles)

data = get(handles.table_objectProperties,'data');
rows = data(:,1);
indx = eventdata.Indices(1);
val  = eventdata.NewData;
virmenEventHandler('changeObjectProperties',{rows{indx},val});

function objectPropertiesSelect_Callback(hObject, eventdata, handles)

highlightFigure;


function objectLocationsEdit_Callback(hObject, eventdata, handles)

data = get(handles.table_objectLocations,'data');
virmenEventHandler('changeObjectLocations',data(:,2:end));

function objectLocationsSelect_Callback(hObject, eventdata, handles)

set(hObject,'userdata',eventdata.Indices(:,1));
highlightFigure;


% --- Executes on selection change in pop_shape.
function pop_shape_Callback(hObject, eventdata, handles)

val = get(hObject,'value');
virmenEventHandler('selectShape',{val,'none',1});


% --- Executes during object creation, after setting all properties.
function pop_shape_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_setShapeColor.
function push_setShapeColor_Callback(hObject, eventdata, handles)

virmenEventHandler('selectShape','n/a')


% --- Executes on button press in push_shapeSymbolic.
function push_shapeSymbolic_Callback(hObject, eventdata, handles)

virmenEventHandler('changeSymbolicShapeLocations',[]);


% --- Executes on button press in push_renameShape.
function push_renameShape_Callback(hObject, eventdata, handles)

virmenEventHandler('renameShape',[]);


function shapePropertiesEdit_Callback(hObject, eventdata, handles)

data = get(handles.table_shapeProperties,'data');
rows = data(:,1);
indx = eventdata.Indices(1);
val  = eventdata.NewData;
virmenEventHandler('changeShapeProperties',{rows(indx),{val}});


function shapePropertiesSelect_Callback(hObject, eventdata, handles)

highlightFigure;



function shapeLocationsEdit_Callback(hObject, eventdata, handles)

data = get(handles.table_shapeLocations,'data');
virmenEventHandler('changeShapeLocations',data(:,2:end));


function shapeLocationsSelect_Callback(hObject, eventdata, handles)

set(hObject,'userdata',eventdata.Indices(:,1));
highlightFigure;


function variablesEdit_Callback(hObject, eventdata, handles)

data = get(handles.table_variables,'data');
rows = data(:,1);
indx = eventdata.Indices(1);
val  = eventdata.NewData;
virmenEventHandler('changeVariables',{rows{indx},val});

function variablesSelect_Callback(hObject, eventdata, handles)

set(hObject,'userdata',eventdata.Indices(:,1));
highlightFigure;


% --- Executes on selection change in pop_experimentCode.
function pop_experimentCode_Callback(hObject, eventdata, handles)

str = get(hObject,'string');
val = get(hObject,'value');
if strcmp(str{val}(1),'[')
    fn = [];
else
    fn = str2func(str{val});
end
    
virmenEventHandler('changeExperimentProperties',{'experimentCode',fn});


% --- Executes during object creation, after setting all properties.
function pop_experimentCode_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_editTransformation.
function push_editTransformation_Callback(hObject, eventdata, handles)

virmenEventHandler('editCode','transformationFunction');


% --- Executes on button press in push_windowManager.
function push_windowManager_Callback(hObject, eventdata, handles)

virmenEventHandler('changeExperimentProperties','n/a');


% --- Executes on button press in push_editMovement.
function push_editMovement_Callback(hObject, eventdata, handles)

virmenEventHandler('editCode','movementFunction');


% --- Executes on button press in push_sortShapes.
function push_sortShapes_Callback(hObject, eventdata, handles)

virmenEventHandler('sortShapes','n/a');

% --- Executes on button press in push_sortObjects.
function push_sortObjects_Callback(hObject, eventdata, handles)

virmenEventHandler('sortObjects','n/a');
