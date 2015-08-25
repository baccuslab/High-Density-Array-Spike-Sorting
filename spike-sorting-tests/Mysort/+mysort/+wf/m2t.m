
function T = m2t(M, nC)
    % converts the multi channel waveforms stores as rows in M (with
    % multiple channels in multiple rows) into the tensor representation.
    % If nS is the number of waveforms, Tf is the number of samples per
    % individual channel and nC is the number of channels, 
    %                 M is (nS*nC) x Tf
    % and the resulting 
    %                 T is Tf x nC x nS
    [nSnC, Tf] = size(M);
    nS = nSnC/nC;
    assert(round(nS)==nS, 'Number of row of M does not match nC!');
    
    T = reshape(M', Tf, nC, nS);