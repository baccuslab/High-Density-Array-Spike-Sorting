function varargout = dataViewerAxes(varargin)
% DATAVIEWERAXES M-file for dataViewerAxes.fig
%      DATAVIEWERAXES, by itself, creates a new DATAVIEWERAXES or raises the existing
%      singleton*.
%
%      H = DATAVIEWERAXES returns the handle to a new DATAVIEWERAXES or the handle to
%      the existing singleton*.
%
%      DATAVIEWERAXES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DATAVIEWERAXES.M with the given input arguments.
%
%      DATAVIEWERAXES('Property','Value',...) creates a new DATAVIEWERAXES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dataViewerAxes_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dataViewerAxes_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dataViewerAxes

% Last Modified by GUIDE v2.5 01-Apr-2011 19:14:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dataViewerAxes_OpeningFcn, ...
                   'gui_OutputFcn',  @dataViewerAxes_OutputFcn, ...
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


% --- Executes just before dataViewerAxes is made visible.
function dataViewerAxes_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dataViewerAxes (see VARARGIN)

% Choose default command line output for dataViewerAxes
handles.output = hObject;

callbacks = {'slider', 'zoomWindowCheckbox', 'centralMsButton', 'zoomPopupmenu'};
for i=1:length(callbacks)    
    set(handles.(callbacks{i}), 'Callback', @(x,y) figutil.evoke_callback('dataaxes', callbacks{i}, x, y));
end
set(handles.dataAxes, 'ButtonDownFcn', @(x,y) figutil.evoke_callback('dataaxes', 'dataAxes', x, y));
set(handles.figure1, 'ButtonDownFcn',  @(x,y) figutil.evoke_callback('dataaxes', 'figure', x, y));
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dataViewerAxes wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dataViewerAxes_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider_Callback(hObject, eventdata, handles)
% hObject    handle to slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
callBackFunction = get(hObject, 'UserData');
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in zoomWindowCheckbox.
function zoomWindowCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to zoomWindowCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of zoomWindowCheckbox
callBackFunction = get(hObject, 'UserData');

% --- Executes on button press in centralMsButton.
function centralMsButton_Callback(hObject, eventdata, handles)
% hObject    handle to centralMsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
callBackFunction = get(hObject, 'UserData');


function centralMsEdit_Callback(hObject, eventdata, handles)
% hObject    handle to centralMsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of centralMsEdit as text
%        str2double(get(hObject,'String')) returns contents of centralMsEdit as a double
callBackFunction = get(hObject, 'UserData');

% --- Executes during object creation, after setting all properties.
function centralMsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to centralMsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in zoomPopupmenu.
function zoomPopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to zoomPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns zoomPopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from zoomPopupmenu
callBackFunction = get(hObject, 'UserData');

% --- Executes during object creation, after setting all properties.
function zoomPopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zoomPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over axes background.
function dataAxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to dataAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over figure background.
function figure1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


