
function sortingErrors(R, sorter, errorType,  varargin)
    P.cutleft = 100;
    P.cutright = 100;
    P.unitList = [];
    P.maxErrorsToPlot = 5;
    P = mysort.util.parseInputs(P, 'sortingErrors', varargin);
    
    
    errorSt = mysort.spiketrain.getErrorSpikeTrain(R, errorType);
    if isempty(P.unitList)
        P.unitList = 1:length(errorSt);
    end
    for i=1:length(errorSt)
        if ~any(P.unitList == i);
            continue
        end
        samples = errorSt{i};
        
        if length(samples)>P.maxErrorsToPlot
            warning('Too many errors to plot (%d), plotting only first %d!',length(samples),P.maxErrorsToPlot);
        end
        for k=1:min(P.maxErrorsToPlot,length(samples))
            s1 = max(1,samples(k)-P.cutleft);
            s2 = samples(k)+P.cutright+sorter.Tf;
            sorter.plotSorting('start', s1, 'stopp', s2, 'eval', R);
            mysort.plot.figureTitle(sprintf('Errortype: %s Unit: %d Error %d (of %d)', errorType, i, k, length(samples)));
            mysort.plot.figureName('sorting, errors');
        end
    end
end