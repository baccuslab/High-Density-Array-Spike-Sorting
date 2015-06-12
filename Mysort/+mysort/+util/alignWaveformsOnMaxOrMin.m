
function [tau X] = alignWaveformsOnMaxOrMin(X, nC, varargin)
warning('This function is depricated! Use mysort.wf.* instead!');
    P.debug = 0;
    P.restrictToIdx = [];
    P.maxIdx = [];
    P.truncate = 1;
    P.axesHandles = [];
    P = mysort.util.parseInputs(P, 'alignWaveformsOnMaxOrMin', varargin);
    
    Tf = size(X,2)/nC;
    assert(round(Tf) == Tf, 'Channel number mismatch!');
    
    if isempty(P.restrictToIdx)
        P.restrictToIdx = 1:Tf;
    end
    
    if ~isempty(P.axesHandles)
        assert(length(P.axesHandles) == 4, 'Need 4 axes handles for this plot!');
        P.debug = 1;
    end
    
    if P.debug && isempty(P.axesHandles);
        figure('color','w');
        P.axesHandles(1) = subplot(2,3,1:2);
        P.axesHandles(2) = subplot(2,3,3);
        P.axesHandles(3) = subplot(2,3,3:4);
        P.axesHandles(4) = subplot(2,3,5);
    end
    
    alignthis = mysort.wf.v2m(X, nC);
    setZeroIdx = setdiff(1:Tf, P.restrictToIdx);
    alignthis(:, setZeroIdx) = 0;
    alignthis = mysort.wf.m2v(alignthis, nC);
    if abs(min(alignthis(:))) > max(alignthis(:))
        alignthis = -alignthis;
    end
    
    if P.debug
        plot(P.axesHandles(1), X'); 
        plot(P.axesHandles(2), alignthis'); 
    end    

    [tau alignthis] = mysort.util.alignWaveformsOnMax(alignthis, nC,...
                        'truncate', P.truncate, 'maxIdx', P.maxIdx);

    X = mysort.util.shiftMCRows(X, tau, nC, 1);

    if P.debug
        plot(P.axesHandles(3), X'); 
        plot(P.axesHandles(4), alignthis'); 
    end       
    