function consoleOutput(handles, str, color)
    old_str = get(handles.Console, 'String');
    if ~iscell(old_str)
        old_str = {old_str};
    end
    maxL = 2000;
    if length(old_str) > maxL
        old_str = old_str(end-maxL+1:end);
    end
    new_string = [old_str; str];
    set(handles.Console, 'ForegroundColor', color);
    set(handles.Console, 'String', new_string);
    set(handles.Console, 'Value', length(new_string));
    drawnow 