
function [tau S] = alignWaveformsInterpolatedMin(S,nC, varargin)
warning('This function is depricated! Use mysort.wf.* instead!');
    if nargin < 2
        nC = 1;
    end
    [tau S] = mysort.util.alignWaveformsInterpolatedMax(-S, nC, varargin{:});
    S = -S;
    