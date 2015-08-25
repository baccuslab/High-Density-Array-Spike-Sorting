function T = v2t(V, nC)
    % converts the multichannel waveforms stored as rows in v with single
    % channels concatenated into the tensor representation
    
    [nS TfnC nSubsampleTaus] = size(V);
    Tf = TfnC/nC;
    assert(round(Tf)==Tf, 'nC does not match dimensions of V!');
    P = permute(V, [2 1 3]);
    T = reshape(P, Tf, nC, nS, nSubsampleTaus);