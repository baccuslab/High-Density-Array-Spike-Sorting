
function [tau A] = alignWaveformsUpsampleMaxOrMin(S, nC, varargin)
warning('This function is depricated! Use mysort.wf.* instead!');
    P.maxIdx = [];
    P.UP = 3;
    P.downsample = 1;
    P = mysort.util.parseInputs(P, 'alignWaveformsUpsampleMaxOrMin', varargin);

    Sup = mysort.util.resampleTensor(mysort.util.m2t(S, nC), P.UP,1);
    Sup = mysort.util.t2m(Sup);
    [tau Sup] = mysort.util.alignWaveformsOnMaxOrMin(Sup, nC, 'maxIdx', P.UP*P.maxIdx);
    if P.downsample
        tau = tau/P.UP;
        A = mysort.util.resampleTensor(mysort.util.m2t(Sup, nC), 1, P.UP);
        A = mysort.util.t2m(A);
    else
        A = Sup;
    end