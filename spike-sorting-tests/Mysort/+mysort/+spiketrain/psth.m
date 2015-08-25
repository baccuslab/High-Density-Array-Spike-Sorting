function [C edges] = psth(spiketrains, binWidth, edges)
    % Computes the psth for all spike trains in cell array spiketrains.
    % Each column of spiketrains is treated individually
    [nT nC] = size(spiketrains);
    if nargin == 2 || isempty(edges);
        st = spiketrains(:);
        st(cellfun(@isempty, st)) = [];
        minTimes = cellfun(@(x) min(x), st);
        maxTimes = cellfun(@(x) max(x), st);
        totalMin = min(minTimes(:));
        totalMax = max(maxTimes(:));
        totalTime = totalMax - totalMin;
        nBins = ceil(totalTime/binWidth);
        edges = totalMin + (0:nBins)*binWidth;
    end
    nEdges = length(edges);
    C = zeros(nEdges, nC);
    
    for c=1:nC
        C(:,c) = histc(cell2mat(spiketrains(:,c)'), edges);
    end