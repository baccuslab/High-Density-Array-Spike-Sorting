
function P = clustering(X, classes, templates, template_ids, varargin)
    P.errorTypes = [];
    P.width = 930;
    P.height = 300;
    P.axesHandles = [];
    P.markerSize = 8;
    P.fh = [];
    P = mysort.util.parseInputs(P, varargin, 'error');
    if isempty(classes)
        classes = ones(size(X,1),1);
    end
    if ndims(classes) == 2
        classes = classes(:,1);
    end
    if ~exist('template_ids', 'var') || isempty(template_ids)
        template_ids = unique(classes);
    end
    if exist('templates', 'var') && ~isempty(templates)
        if nargin < 4 || isempty(template_ids)
            template_ids = 1: size(templates,1);
        else
            assert(length(template_ids) == size(templates,1), 'There must be as many template IDs as templates!');
        end
    else
        templates = mysort.util.calculateClassMeans(X, classes);
    end

    def = mysort.util.defs();
    if isempty(P.errorTypes)
        P.errorTypes = ones(size(X,1), 1) * def.tp;
    end
    
    nClasses = length(template_ids);
    dims = size(X,2);
    if dims == 3 || dims == 6
        nPlots = dims/3;
        plot3d = 1;
    else
        plot3d = 0;
        nPlots = ceil(dims/2);
    end
     
    if isempty(P.axesHandles)
        if isempty(P.fh)
            P.fh = mysort.plot.figure('width', P.width, 'height', P.height);
        end
        mysort.plot.figureName('Clustering');
        ax = mysort.plot.subplots([1 nPlots]);
    else
        ax = P.axesHandles;
    end
    cW = 5;
    pW = P.markerSize;
    cLW = 2;
    for p=1:nPlots
        axes(ax(p)); hold on
        for i=1:nClasses
            myId = template_ids(i);
            idx   = classes==myId;
            idxCL = classes==myId & P.errorTypes~=def.tp & P.errorTypes~=def.fp;
            idxFP = classes==myId & P.errorTypes==def.fp;
            if plot3d
                plot3(X(idx,3*p-2), X(idx,3*p-1), X(idx,3*p), '.', 'color', mysort.plot.vectorColor(myId), 'markersize', pW);
                plot3(X(idxCL,3*p-2), X(idxCL,3*p-1), X(idxCL,3*p), 'o', 'color', [1 0 0], 'markersize', cW, 'linewidth', cLW);
                plot3(X(idxFP,3*p-2), X(idxFP,3*p-1), X(idxFP,3*p), 'o', 'color', [0 0 1], 'markersize', cW, 'linewidth', cLW);
            else
                plot(X(idx,2*p-1), X(idx,2*p), '.', 'color', mysort.plot.vectorColor(myId), 'markersize', pW);
                plot(X(idxCL,2*p-1), X(idxCL,2*p), 'o', 'color', [1 0 0], 'markersize', cW, 'linewidth', cLW);
                plot(X(idxFP,2*p-1), X(idxFP,2*p), 'o', 'color', [0 0 1], 'markersize', cW, 'linewidth', cLW);
            end
            
            if plot3d
                plot3(templates(i,3*p-2), templates(i,3*p-1), templates(i,3*p), 'x', 'color', 'k', 'markersize', 14, 'linewidth', 4);
                plot3(templates(i,3*p-2), templates(i,3*p-1), templates(i,3*p), 'x', 'color', mysort.plot.vectorColor(myId), 'markersize', 14, 'linewidth', 2);
            else
                plot(templates(i,2*p-1), templates(i,2*p), 'x', 'color', 'k', 'markersize', 14, 'linewidth', 4);
                plot(templates(i,2*p-1), templates(i,2*p), 'x', 'color', mysort.plot.vectorColor(myId), 'markersize', 14, 'linewidth', 2);
%                 plot(templates(i, 2*p-1), templates(i, 2*p), 'x', 'color', [0 0 0], 'markersize', 14, 'linewidth', 4);
%                 plot(templates(i, 2*p-1), templates(i, 2*p), 'x', 'color', mysort.plot.vectorColor(myId), 'markersize', 12, 'linewidth', 2);
            end
        end
%         set(ax(p), 'color', 'k');
    end
    P.axesHandles = ax;