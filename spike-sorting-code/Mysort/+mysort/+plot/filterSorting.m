
function filterSorting(X, gdf, T, D, varargin)
    P.figureHandle = [];
    P.axesHandle = [];
    P.sampleOffset = [];
    P.srate = [];
    P.Y = [];
    P.spacer = [];
    P = mysort.util.parseInputs(P, 'filterSorting', varargin);
    
    if isempty(P.axesHandle)
        if isempty(P.figHandle)
            P.figHandle = mysort.plot.figure('color','w');
        end
        P.axesHandle = mysort.plot.subplots(2,1);
    end
    mysort.plot.figureName('FilterSorting'); 

    mysort.plot.sorting(X,gdf,T,'axesHandle', P.axesHandle(1),...
                                'spacer', P.spacer, 'srate', P.srate,...
                                'sampleOffset', P.sampleOffset);
    