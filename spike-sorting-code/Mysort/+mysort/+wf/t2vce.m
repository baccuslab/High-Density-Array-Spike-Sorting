function vce = t2vce(t)
    [Tf, nC, nT] = size(t);
    
    if isempty(t)
        vce = [];
        return
    end
    
    vce = mysort.wf.vte2vce(mysort.wf.t2v(t), nC);