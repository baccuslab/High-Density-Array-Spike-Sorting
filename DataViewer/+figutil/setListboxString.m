function setListboxString(lb, str)
    if iscell(str)
        L = length(str);
    else
        L = 1;
    end
    val = get(lb, 'Value');
    if any(val > L) || any(val<1)
        newval = val(val<=L & val>=1);
        if isempty(newval)
            newval = 1;
        end
        set(lb, 'Value', newval);
    end
    set(lb, 'String', str);

    