function varargout = hdsort(varargin)
% HDSORT MATLAB code for hdsort.fig
%      HDSORT, by itself, creates a new HDSORT or raises the existing
%      singleton*.
%
%      H = HDSORT returns the handle to a new HDSORT or the handle to
%      the existing singleton*.
%
%      HDSORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HDSORT.M with the given input arguments.
%
%      HDSORT('Property','Value',...) creates a new HDSORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before hdsort_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to hdsort_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help hdsort

% Last Modified by GUIDE v2.5 13-Mar-2012 19:03:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @hdsort_OpeningFcn, ...
                   'gui_OutputFcn',  @hdsort_OutputFcn, ...
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


% --- Executes just before hdsort is made visible.
function hdsort_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to hdsort (see VARARGIN)

% Choose default command line output for hdsort

handles.output = hObject;
handles = hdmeagui.init(handles);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes hdsort wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = hdsort_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in UpdateButton.
function UpdateButton_Callback(hObject, eventdata, handles)
% hObject    handle to UpdateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function ProjThreshEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ProjThreshEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ProjThreshEdit as text
%        str2double(get(hObject,'String')) returns contents of ProjThreshEdit as a double


% --- Executes during object creation, after setting all properties.
function ProjThreshEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ProjThreshEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in TemplateListbox.
function TemplateListbox_Callback(hObject, eventdata, handles)
% hObject    handle to TemplateListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns TemplateListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TemplateListbox


% --- Executes during object creation, after setting all properties.
function TemplateListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TemplateListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PlotTemplateButton.
function PlotTemplateButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotTemplateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function CutSpikesTf_Callback(hObject, eventdata, handles)
% hObject    handle to CutSpikesTf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CutSpikesTf as text
%        str2double(get(hObject,'String')) returns contents of CutSpikesTf as a double


% --- Executes during object creation, after setting all properties.
function CutSpikesTf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CutSpikesTf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CutSpikesCutleft_Callback(hObject, eventdata, handles)
% hObject    handle to CutSpikesCutleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CutSpikesCutleft as text
%        str2double(get(hObject,'String')) returns contents of CutSpikesCutleft as a double


% --- Executes during object creation, after setting all properties.
function CutSpikesCutleft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CutSpikesCutleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DisplaySpikesTf_Callback(hObject, eventdata, handles)
% hObject    handle to DisplaySpikesTf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DisplaySpikesTf as text
%        str2double(get(hObject,'String')) returns contents of DisplaySpikesTf as a double


% --- Executes during object creation, after setting all properties.
function DisplaySpikesTf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DisplaySpikesTf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DisplaySpikesCutleft_Callback(hObject, eventdata, handles)
% hObject    handle to DisplaySpikesCutleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DisplaySpikesCutleft as text
%        str2double(get(hObject,'String')) returns contents of DisplaySpikesCutleft as a double


% --- Executes during object creation, after setting all properties.
function DisplaySpikesCutleft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DisplaySpikesCutleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SpikeDetectionThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to SpikeDetectionThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SpikeDetectionThreshold as text
%        str2double(get(hObject,'String')) returns contents of SpikeDetectionThreshold as a double


% --- Executes during object creation, after setting all properties.
function SpikeDetectionThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SpikeDetectionThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NPlotChans_Callback(hObject, eventdata, handles)
% hObject    handle to NPlotChans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NPlotChans as text
%        str2double(get(hObject,'String')) returns contents of NPlotChans as a double


% --- Executes during object creation, after setting all properties.
function NPlotChans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NPlotChans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AcceptTemplateButton.
function AcceptTemplateButton_Callback(hObject, eventdata, handles)
% hObject    handle to AcceptTemplateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in SaveButton.
function SaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in NewTemplateButton.
function NewTemplateButton_Callback(hObject, eventdata, handles)
% hObject    handle to NewTemplateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in UpsampleButton.
function UpsampleButton_Callback(hObject, eventdata, handles)
% hObject    handle to UpsampleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in AlignButton.
function AlignButton_Callback(hObject, eventdata, handles)
% hObject    handle to AlignButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
