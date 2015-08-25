
function Y = shiftMCRows(X, tau, nC, trunc)
    % shifts the rows of X by tau. X can contain multichannel data, with
    % all rows concatinated. If so, the number of channels has to be given
    % by nC
    [nRows inDim] = size(X);
    if nargin < 4; trunc = 0; end
    if nargin < 3; nC = 1; end
    dimPerChan = inDim/nC;
    
    assert(round(dimPerChan) == dimPerChan, 'Number of channel does not match X!');
    assert(trunc==0 || trunc == 1, 'trunc has to be true or false!');
    if size(tau,2)>1; tau = tau'; end
    assert(size(tau,2) == 1, 'tau must be a vector!');
    
    if nC == 1
        Y = mysort.util.shiftRows(X, tau, trunc);
        return;
    end

    XM = mysort.wf.v2m(X, nC);
    tauM = repmat(tau',nC,1);
    tauM = tauM(:);
    XMS = mysort.util.shiftRows(XM, tauM, trunc);
    Y = mysort.wf.m2v(XMS, nC);
    