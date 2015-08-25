function T = tShift(T, tau, trunc)
    % shifts the single channel waveforms stored in T
    % truncates in a way, that A has the same dimensions as T
    % If a waveform has a shift of 0, and trunc is 1 it will be unchanged.
    
    if nargin < 3; trunc = 0; end
    tau = tau(:);
    assert(~any(tau~=round(tau)), 'tau must contain only integers!');
    nC = size(T,2);
    tau = repmat(tau', nC, 1);
    tau = tau(:);    
    T = mysort.wf.t2m(T);
    T = mysort.util.shiftRows(T, tau, trunc);
    T = mysort.wf.m2t(T, nC);   