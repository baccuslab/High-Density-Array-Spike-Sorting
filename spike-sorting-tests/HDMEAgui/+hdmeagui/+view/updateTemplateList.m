function updateTemplateList(handles)
    h = handles.TemplateListbox;    
    T = handles.GUI.T;
    tid = handles.GUI.selTemplateId;
    tidx = T.tId2tIdx(tid);
    str = {};
    I = T.getIterator();
    while I.hasNext()
        t = I.next();
        name = 'T';
        if ~t.isAccepted
            name = 't';
        end
        str{I.idx} = sprintf('%s %d (# %5d)', name, t.id, t.getNWfs());
    end
    set(h, 'string', str, 'value', tidx);
   
    