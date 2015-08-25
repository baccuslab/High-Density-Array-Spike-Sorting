function handles = ChannelAmpAxes(hObject, eventdata, handles)
    point = get(gca, 'CurrentPoint');
    chan = round(point(1,1));
    thr = point(1,2);
    handles.GUI.setChanAmpThr(chan, thr);
    if isfield(handles, 'channelAmpHandle')
        try
            delete(handles.channelAmpHandle);
        catch
        end
    end
    set(handles.ChannelAmpAxes, 'nextplot', 'add');
    handles.channelAmpHandle = ...
        plot(chan+[-.5 .5], thr+[0 0], ...
        '-r', 'linewidth', 2);