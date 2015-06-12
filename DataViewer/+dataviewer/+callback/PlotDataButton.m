function handles = PlotDataButton(hObject, eventdata, handles)
    P = dataviewer.util.getHandleValues(handles);
    
    % handles.Axes = dataViewerAxes();
    if length(P.trialIDs) ~= 1; handles.warning_func('Only a single trial can be plotted at the moment!'); return; end
    % Ignore selected events and plot all
    P.eventIDs = [];
    P.events = [];
    E = handles.DH.getEventsStruct('otherP', P);
    % remove lastSt Event
    for i=length(E):-1:1
        if strcmp(E(i).name, 'lastSt')
            E(i) = [];
        end
    end
    [X ] = handles.DH.getData('otherP', P); %, 'downsampleTo', 3000
    srate = handles.DH.getSamplesPerSecond();
    % handles.dataAxesSpacer = mysort.plot.mc(X,'axesHandle',handles.dataAxes,'figure',0, 'srate', 1/32, 'color', {'k'},...
    %     'events',E,'eventsfirst',1);
    if ~isempty(P.from); offset = P.from; else offset = 0; end
    handles.dataAxesSpacer = mysort.plot.mc(X, 'srate', 1000/srate, 'color', {'k'},...
         'events',E,'eventsfirst', 1, 'sampleOffset', offset); %, 'xticks', T/32);
    handles.dataAxes = gca;
    handles.new_figure_handles = gcf;
    xlabel('time [ms]');
    mysort.plot.figureTitle(P.figureTitle);
    