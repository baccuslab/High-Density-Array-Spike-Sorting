function V = vSubsel(V, nC, idx)
    % subselects the time index into multi channel waveforms on every
    % channel if the waveforms are stored in concatenated form as rows of
    % matrix V
    if isempty(V)
        return
    end
    vidx = mysort.wf.vSubIdx(V,nC,idx);
    V = V(:,vidx);