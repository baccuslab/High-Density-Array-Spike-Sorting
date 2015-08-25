function Vs = vShift(V, nC, tau, trunc)
    % shifts the single channel waveforms stored in a concatenated way in
    % V. Ever row of V is multi channel waveforms
    % truncates the sides of the resulting shifted waveforms if trunc is
    % 1. truncates in a way, that V and Vs have the same dimensions.
    % If a waveform has a shift of 0, and trunc is 1 it will be unchanged.
    
    assert(nargin >= 3, 'V and tau need to be provided!');
    if nargin < 4; trunc = 0; end
    [nS, nCTf] = size(V);
    assert(nS == length(tau), 'For every row of V a shift value needs to be provided!');
    assert(trunc==0 || trunc==1, 'Trunc must be either true or false!');
    assert(size(tau,1) == 1 || size(tau,2) == 1, 'tau must not be a matrix!');
    if size(tau,1)>size(tau,2); tau = tau'; end
    assert(~any(tau~=round(tau)), 'tau must contain only integers!');
    Tf = nCTf/nC;
    assert(round(Tf) == Tf, 'nC does not match dimensions of V!');
    
    tau_ = repmat(tau, nC, 1);
    tau_ = tau_(:);
    M = mysort.wf.v2m(V,nC);
    Ms = mysort.util.shiftRows(M, tau_, trunc);
    Vs = mysort.wf.m2v(Ms, nC);   