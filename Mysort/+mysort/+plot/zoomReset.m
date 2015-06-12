function zoomReset(axesHandle)
    try
        setappdata(axesHandle, 'matlab_graphics_resetplotview', []);
    catch
    end