function handles = upsampleSpikes(handles)
    D = handles.DATA;
    t = handles.INTERACTIONS.tIdx;
    if isempty(D.templates.cutSpikes(t))
        disp('no spikes');
        return
    end
    sp = D.templates.cutSpikes(t);
    for i=1:size(sp,1)
        
    end