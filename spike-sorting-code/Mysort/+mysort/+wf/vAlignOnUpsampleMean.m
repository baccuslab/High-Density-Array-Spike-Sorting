function [A, tau] = vAlignOnUpsampleMean(S, nC, varargin)
    A = mysort.wf.v2t(S,nC);
    [A, tau] = mysort.wf.tAlignOnUpsampleMean(A, varargin{:});
    A = mysort.wf.t2v(A);