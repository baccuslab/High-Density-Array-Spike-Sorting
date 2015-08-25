
function M = v2m(V, nC)
    % converts the vector represention V into the matrix representation M.
    % single channel are concatenated in V so that every row of V
    % represents one multichannel waveform. In M the single channel
    % waveforms are below each other
  
    [nS nCTf] = size(V);
    Tf = nCTf/nC;
    assert(round(Tf) == Tf, 'number of columns in V does not match nC!');
    
    M = reshape(V', Tf, nS*nC)';
end