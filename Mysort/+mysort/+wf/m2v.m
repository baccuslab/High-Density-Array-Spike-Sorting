function V = m2v(M, nC)
    % converts waveforms in matrix representation (i.e. every row of M is a
    % single channel waveform. if a waveform has multiple channels these
    % are in individual rows below each other in m) into the vector
    % representation where the single channel waveforms are concatenated
    %
    % if nC is not supplied it is assumed, that m contains only one
    % waveform and all rows of m are concatenated
    if nargin == 1
        % opposite of v2m: x := v2m(m2v(x),size(x,1))
        V = M';
        V = V(:)';
        return
    end
    
    nS = size(M,1)/nC;
    assert(round(nS) == nS, 'Number of rows of M must be a multiple of nC!');
    Tf = size(M,2);

    V = reshape(M', Tf*nC, nS)';