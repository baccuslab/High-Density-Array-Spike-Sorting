function [C P] = getNoiseCovariance(X, varargin)
    P = mysort.configuration();
    P = mysort.util.parseInputs(P, varargin, 'merge');
    
    X = mysort.checkData(X, P.srate);