
function count = figureChildrenSet(handle, parameter, value)
    count = 0;
    try
        set(handle, parameter, value);
        count = count +1;
    catch
        % Ignore this
    end
    c = get(handle, 'Children');
    for i=1:length(c)
%         list = get(c(i));
%         if isfield(list, 'Tag') && ~strcmp(list.Tag, 'legend');
%             %ignore
%         else
            count = count + mysort.plot.figureChildrenSet(c(i), parameter, value);
%         end
    end