% TODO-
%{
-update slider position when edit text change
-don't let edit1 exceed the num_frames
-Show text box from very start
- fix hard values for contrast, and median
-filter_flag needs to be in handles for all to access,
-fitler_flag on slider change
%}



function varargout = pupil_gui(varargin)
% PUPIL_GUI MATLAB code for pupil_gui.fig
%      PUPIL_GUI, by itself, creates a new PUPIL_GUI or raises the existing
%      singleton*.
%
%      H = PUPIL_GUI returns the handle to a new PUPIL_GUI or the handle to
%      the existing singleton*.
%
%      PUPIL_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PUPIL_GUI.M with the given input arguments.
%
%      PUPIL_GUI('Property','Value',...) creates a new PUPIL_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pupil_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pupil_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pupil_gui

% Last Modified by GUIDE v2.5 30-Jan-2017 10:28:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pupil_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @pupil_gui_OutputFcn, ...
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




% --- Executes just before pupil_gui is made visible.
function pupil_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pupil_gui (see VARARGIN)

% Choose default command line output for pupil_gui
handles.output = hObject;
handles.cam0=varargin{1};

num_frames=length(handles.cam0);
set(handles.slider1,'value',1, 'min',1, 'max',num_frames,'SliderStep',[1, 1] /(num_frames -1));

handles.filter_flag='NORMAL';

imshow( handles.cam0{1}, 'Parent', handles.axes1 );

% Update handles structure
guidata(hObject, handles);





% UIWAIT makes pupil_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pupil_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
%slider_value = int32(get(, 'Value')); gets slider value
%anywher
slider_value=round(get(hObject,'Value'));
set(hObject, 'Value', slider_value);

% img=handles.cam0{slider_value};
% cla;
% imshow(img , 'Parent', handles.axes1 );
% set edit1
set(handles.edit1, 'String', num2str(slider_value));

apply_image_filter();


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit1_value=str2double(get(hObject,'String'));
set(handles.slider1,'value',edit1_value)
apply_image_filter();
% img=handles.cam0{edit1_value};
% cla;
% imshow(img , 'Parent', handles.axes1 );


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in radiobutton7.
function radiobutton7_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.filter_flag='CONTRAST';
guidata(gcbo,handles);
%show image, just once for radio button push
apply_image_filter();


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)

handles.filter_flag='MEDIAN';
guidata(gcbo,handles);
apply_image_filter();
%show image, just once for radio button push



function img = apply_image_filter()

handles = guidata(gcbo);

filter_flag=handles.filter_flag;
cam0=handles.cam0{str2num(get(handles.edit1,'String'))};
    switch filter_flag
        case 'NORMAL'
            img=cam0;
        case 'MEDIAN'
            img=medfilt2(cam0,[7 7]);
        case 'CONTRAST'
            range_for_min={50:183,20:155}; %range to look for pixels values, (1)=x, (2)=y;  250x250 = {50:200,50:200} 
            img1=medfilt2(cam0,[7 7]);
            img2=img1;
            img2(find(img1<=min(min(img1(range_for_min{1},range_for_min{2})))+13))=0;
            img=img2;
            
    end
    cla; imshow(img , 'Parent', handles.axes1 );
    
    guidata(gcbo,handles); % update the handles


    


% --- Executes on button press in radiobutton8.
function radiobutton8_Callback(hObject, eventdata, handles)
handles.filter_flag='NORMAL';
guidata(gcbo,handles);
%show image, just once for radio button push
apply_image_filter();

