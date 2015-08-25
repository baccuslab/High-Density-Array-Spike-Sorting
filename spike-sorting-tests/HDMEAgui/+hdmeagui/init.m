function handles = init(handles)
    % DEFINE constants
    CONFIG = struct('color_selected_spikes', 'g', ...
                    'color_unselected_spikes', [.4 .4 .4], ...
                    'color_selected_template', [.9 .5 .0], ...
                    'color_unselected_template', 'k', ...
                    'color_normhist_new_thresh', 'c', ...
                    'color_normhist_old_thresh', 'k');
    % define numeric values settable by the GUI
    GUI_CONFIG = struct('CutSpikesTf', 30, ...
                        'CutSpikesCutleft', 9, ...
                        'DisplaySpikesTf', 15, ...
                        'DisplaySpikesCutleft', 5,...
                        'NPlotChans', 10);            
    % get user data from gui call
    UD = get(handles.figure1, 'UserData');
    handles = hdmeagui.initInput(handles, UD, CONFIG, GUI_CONFIG);
    
    handles.GUI = hdmeagui.GuiManager(handles.WF, handles.T);

    % check if GUI elements for parameter values exist and set init values
    fnames = fieldnames(handles.GUI_CONFIG);
    postCallbackCheck = @hdmeagui.postCallbackCheck;
    for i=1:length(fnames)
        str = fnames{i};
        val = handles.GUI_CONFIG.(str);    
        assert(isfield(handles, str), ...
            sprintf('GUI parameter (%s) could not be associated with GUI element!', str));
        if isscalar(handles.GUI_CONFIG.(fnames{i}))
            cback = @(x,y) figutil.evoke_callback_editFieldSync(str, @str2num, x, y, postCallbackCheck);
            set(handles.(str), 'string', num2str(val), 'callback', cback);
        end
    end


    % define and register callbacks
    callbacks = {'UpdateButton', 'AcceptTemplateButton', 'SaveButton',...
                 'PlotTemplateButton', 'NewTemplateButton', 'UpsampleButton', 'AlignButton'};
    postCallbackCheck = @hdmeagui.postCallbackCheck;
    for i=1:length(callbacks)    
        set(handles.(callbacks{i}), 'Callback', @(x,y) figutil.evoke_callback('hdmeagui', callbacks{i}, x, y, postCallbackCheck));
    end
    set(handles.NormHistAxes, 'ButtonDownFcn', @(x,y) figutil.evoke_callback('hdmeagui', 'NormHistAxes', x, y, postCallbackCheck));
    set(handles.ChannelAmpAxes, 'ButtonDownFcn', @(x,y) figutil.evoke_callback('hdmeagui', 'ChannelAmpAxes', x, y, postCallbackCheck));
    set(handles.ChannelAmpAxes, 'NextPlot', 'replacechildren');


    handles.GUI.featuresChanged = true;
    handles = handles.GUI.update(handles);

