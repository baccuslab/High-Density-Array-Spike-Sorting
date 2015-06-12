function varargout = dataViewer(varargin)
% DATAVIEWER M-file for dataViewer.fig

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_OpeningFcn, ...
                   'gui_OutputFcn',  @main_OutputFcn, ...
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


% --- Executes just before main is made visible.
function main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main (see VARARGIN)

% Choose default command line output for main
handles.output = hObject;
handles.warning_func = @(x) dataviewer.util.consoleOutput(handles, ['WARNING !!!: ' x], 'r');
handles.debug_func = @(x) dataviewer.util.consoleOutput(handles, ['DataHandle: ' x], 'k');
handles.DH = dataviewer.DataHandle('dbconfig', db.munk.DBconfig_Select, ...
    'log_function', handles.debug_func);

handles.dataAxes = [];
handles.fh_store = [];
handles.new_figure_handles = [];

callbacks = {'Experiment', 'Block', 'Tetrode', 'Trial', 'Analysis', 'Unit',...
             'Channel', 'SelectTrialsetButton', 'SelectAnalysisIDButton', ...
             'SelectUnitIDButton', 'Trialset', 'trial_ok',...
             'trial_rewarded', 'trial_file_error', 'load_1', 'load_2', ...
             'load_3', 'load_4', 'load_error', 'PrintSelectedParameters', ...
             'CloseAllButton', ...
             'PlotDataButton', 'PlotIntraClusterPCAButton', 'PlotProjectionButton',...
             'PlotClusteringButton','PlotISIButton','PlotXCorrButton',...
             'PlotSortingButton', 'PlotPSTHButton', 'PlotLoadRewardReaction', ...
             'PlotRates', 'PlotSpectogram', ...
             'PlotAmplitudeDriftButton', 'PlotSpikesButton', 'PlotTemplatesButton'};
postCallbackCheck = @dataviewer.util.postCallbackCheck;
for i=1:length(callbacks)    
    set(handles.(callbacks{i}), 'Callback', @(x,y) figutil.evoke_callback('dataviewer', callbacks{i}, x, y, postCallbackCheck));
end
set(handles.AnalysisText, 'ButtonDownFcn', @(x,y) figutil.evoke_callback('dataviewer', 'AnalysisText', x, y));
%set(handles.AnalysisText, 'Callback', @(x,y) figutil.evoke_callback('dataviewer', 'AnalysisText', x, y));

handles.ids.experiments = [];
handles.ids.blocks      = [];
handles.ids.tetrodes    = [];
handles.ids.channels    = [];
handles.ids.trials      = [];
handles.ids.analysis    = [];
handles.ids.units       = [];
handles.selection_names.unit = [];

handles.trialsSelectedBy = 'listbox';

experiments = handles.DH.getExperiments();
figutil.setListboxString(handles.Experiment, experiments(:,2));
handles.ids.experiments = cell2mat(experiments(:,1));

events      = handles.DH.getEventTypes();
figutil.setListboxString(handles.Event, events(:,1));
handles.ids.events = events(:,1);

handles = dataviewer.callback.Experiment(hObject, eventdata, handles);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes main wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in Experiment.
function Experiment_Callback(hObject, eventdata, handles)
% hObject    handle to Experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Experiment contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Experiment


% --- Executes during object creation, after setting all properties.
function Experiment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Block.
function Block_Callback(hObject, eventdata, handles)
% hObject    handle to Block (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Block contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Block



% --- Executes during object creation, after setting all properties.
function Block_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Block (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Trial.
function Trial_Callback(hObject, eventdata, handles)
% hObject    handle to Trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Trial contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Trial

% --- Executes during object creation, after setting all properties.
function Trial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Tetrode.
function Tetrode_Callback(hObject, eventdata, handles)
% hObject    handle to Tetrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Tetrode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Tetrode


% --- Executes during object creation, after setting all properties.
function Tetrode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tetrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Analysis.
function Analysis_Callback(hObject, eventdata, handles)
% hObject    handle to Analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Analysis contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Analysis


% --- Executes during object creation, after setting all properties.
function Analysis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Unit.
function Unit_Callback(hObject, eventdata, handles)
% hObject    handle to Unit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Unit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Unit
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Unit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Unit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in PlotDataButton.
function PlotDataButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in PlotSortingButton.
function PlotSortingButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotSortingButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in PlotSpikesButton.
function PlotSpikesButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotSpikesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in PlotClusteringButton.
function PlotClusteringButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotClusteringButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in PlotIntraClusterPCAButton.
function PlotIntraClusterPCAButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotIntraClusterPCAButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in PlotProjectionButton.
function PlotProjectionButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotProjectionButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in PlotISIButton.
function PlotISIButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotISIButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on selection change in SpikeTrainMode.
function SpikeTrainMode_Callback(hObject, eventdata, handles)
% hObject    handle to SpikeTrainMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns SpikeTrainMode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SpikeTrainMode


% --- Executes during object creation, after setting all properties.
function SpikeTrainMode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SpikeTrainMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in Channel.
function Channel_Callback(hObject, eventdata, handles)
% hObject    handle to Channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Channel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Channel


% --- Executes during object creation, after setting all properties.
function Channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SelectTrialSetButton.
function SelectTrialsetButton_Callback(hObject, eventdata, handles)
% hObject    handle to SelectTrialSetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




function Trialset_Callback(hObject, eventdata, handles)
% hObject    handle to Trialset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Trialset as text
%        str2double(get(hObject,'String')) returns contents of Trialset as a double


% --- Executes during object creation, after setting all properties.
function Trialset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Trialset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PlotDataButton.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to PlotDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in PlotAmplitudeDriftButton.
function PlotAmplitudeDriftButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotAmplitudeDriftButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dataviewer.close(handles);


% --- Executes on button press in PlotXCorrButton.
function PlotXCorrButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotXCorrButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in PlotTemplatesButton.
function PlotTemplatesButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotTemplatesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in SpikeMode.
function SpikeMode_Callback(hObject, eventdata, handles)
% hObject    handle to SpikeMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns SpikeMode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SpikeMode


% --- Executes during object creation, after setting all properties.
function SpikeMode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SpikeMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cutleft_Callback(hObject, eventdata, handles)
% hObject    handle to cutleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cutleft as text
%        str2double(get(hObject,'String')) returns contents of cutleft as a double


% --- Executes during object creation, after setting all properties.
function cutleft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cutleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tf_Callback(hObject, eventdata, handles)
% hObject    handle to Tf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tf as text
%        str2double(get(hObject,'String')) returns contents of Tf as a double


% --- Executes during object creation, after setting all properties.
function Tf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in TemplateModePopup.
function TemplateModePopup_Callback(hObject, eventdata, handles)
% hObject    handle to TemplateModePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns TemplateModePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TemplateModePopup


% --- Executes during object creation, after setting all properties.
function TemplateModePopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TemplateModePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function from_Callback(hObject, eventdata, handles)
% hObject    handle to from (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of from as text
%        str2double(get(hObject,'String')) returns contents of from as a double


% --- Executes during object creation, after setting all properties.
function from_CreateFcn(hObject, eventdata, handles)
% hObject    handle to from (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function to_Callback(hObject, eventdata, handles)
% hObject    handle to to (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of to as text
%        str2double(get(hObject,'String')) returns contents of to as a double


% --- Executes during object creation, after setting all properties.
function to_CreateFcn(hObject, eventdata, handles)
% hObject    handle to to (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function srate_Callback(hObject, eventdata, handles)
% hObject    handle to srate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of srate as text
%        str2double(get(hObject,'String')) returns contents of srate as a double


% --- Executes during object creation, after setting all properties.
function srate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to srate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PlotPSTHButton.
function PlotPSTHButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotPSTHButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in Event.
function Event_Callback(hObject, eventdata, handles)
% hObject    handle to Event (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Event contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Event


% --- Executes during object creation, after setting all properties.
function Event_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Event (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function Console_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Console (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function Binsize_Callback(hObject, eventdata, handles)
% hObject    handle to Binsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Binsize as text
%        str2double(get(hObject,'String')) returns contents of Binsize as a double


% --- Executes during object creation, after setting all properties.
function Binsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Binsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Console.
function Console_Callback(hObject, eventdata, handles)
% hObject    handle to Console (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Console contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Console


% --- Executes on button press in PrintSelectedParameters.
function PrintSelectedParameters_Callback(hObject, eventdata, handles)
% hObject    handle to PrintSelectedParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over AnalysisText.
function AnalysisText_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to AnalysisText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in stacked.
function stacked_Callback(hObject, eventdata, handles)
% hObject    handle to stacked (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of stacked


% --- Executes on button press in trial_ok.
function trial_ok_Callback(hObject, eventdata, handles)
% hObject    handle to trial_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trial_ok


% --- Executes on selection change in trial_rewarded.
function trial_rewarded_Callback(hObject, eventdata, handles)
% hObject    handle to trial_rewarded (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns trial_rewarded contents as cell array
%        contents{get(hObject,'Value')} returns selected item from trial_rewarded


% --- Executes during object creation, after setting all properties.
function trial_rewarded_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trial_rewarded (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in trial_file_error.
function trial_file_error_Callback(hObject, eventdata, handles)
% hObject    handle to trial_file_error (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trial_file_error


% --- Executes on button press in load_1.
function load_1_Callback(hObject, eventdata, handles)
% hObject    handle to load_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of load_1


% --- Executes on button press in load_2.
function load_2_Callback(hObject, eventdata, handles)
% hObject    handle to load_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of load_2


% --- Executes on button press in load_3.
function load_3_Callback(hObject, eventdata, handles)
% hObject    handle to load_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of load_3


% --- Executes on button press in load_4.
function load_4_Callback(hObject, eventdata, handles)
% hObject    handle to load_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of load_4


% --- Executes on button press in load_error.
function load_error_Callback(hObject, eventdata, handles)
% hObject    handle to load_error (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of load_error


% --- Executes on button press in SelectAnalysisIDButton.
function SelectAnalysisIDButton_Callback(hObject, eventdata, handles)
% hObject    handle to SelectAnalysisIDButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function ForcedAnalysisID_Callback(hObject, eventdata, handles)
% hObject    handle to ForcedAnalysisID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ForcedAnalysisID as text
%        str2double(get(hObject,'String')) returns contents of ForcedAnalysisID as a double


% --- Executes during object creation, after setting all properties.
function ForcedAnalysisID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ForcedAnalysisID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SelectUnitIDButton.
function SelectUnitIDButton_Callback(hObject, eventdata, handles)
% hObject    handle to SelectUnitIDButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function ForcedUnitID_Callback(hObject, eventdata, handles)
% hObject    handle to ForcedUnitID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ForcedUnitID as text
%        str2double(get(hObject,'String')) returns contents of ForcedUnitID as a double


% --- Executes during object creation, after setting all properties.
function ForcedUnitID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ForcedUnitID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PlotRates.
function PlotRates_Callback(hObject, eventdata, handles)
% hObject    handle to PlotRates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in PlotLoadRewardReaction.
function PlotLoadRewardReaction_Callback(hObject, eventdata, handles)
% hObject    handle to PlotLoadRewardReaction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in PlotRaster.
function PlotRaster_Callback(hObject, eventdata, handles)
% hObject    handle to PlotRaster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in CloseAllButton.
function CloseAllButton_Callback(hObject, eventdata, handles)
% hObject    handle to CloseAllButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in PlotSpectogram.
function PlotSpectogram_Callback(hObject, eventdata, handles)
% hObject    handle to PlotSpectogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


