function V = vSubChanSel(V, nC, chan_idx)
    % subselects the channel index into multi channel waveforms on every
    % channel specified in chan_idx 
    if numel(V)>5000
        % This is faster for large matrices
        Tf = size(V,2)/nC;
        idx = repmat((1:Tf)', 1, nC) + repmat((0:nC-1), Tf, 1)*Tf;
        idx = idx(:, chan_idx);
        V = V(:, idx(:));
    else
        T = mysort.wf.v2t(V, nC);
        T = T(:,chan_idx,:);
        V = mysort.wf.t2v(T);
    end