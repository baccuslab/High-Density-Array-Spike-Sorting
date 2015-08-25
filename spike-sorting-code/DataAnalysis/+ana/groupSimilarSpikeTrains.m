function [tidx D] = groupSimilarSpikeTrains(gdf, T, meanNoiseStd, varargin)
    error('Dont use this function');
    P.minPercentOfEqualSPikes = 20;
    P.nMaxSpikesPerSpikeTrain = 2000;
    P.binSize = 4;
    P.maxLag = 4;
    P.MES = [];
    [P uP] = mysort.util.parseInputs(P, varargin, 'split');
    duP = mysort.util.deflateP(uP);

    u = unique(gdf(:,1));
    
    %% Group accoring to template similarity
    if nargout > 1
        [groupsT maxT d1] = ana.mergeTemplates(T, meanNoiseStd, duP{:});
    else
        [groupsT maxT] = ana.mergeTemplates(T, meanNoiseStd, duP{:});
    end
    Tmerged1 = T(:,:,maxT);
    gdf(~ismember(gdf(:,1),u(maxT)),:) = [];

    %% Group according to spike train similarity
    groupsSp = num2cell(1:length(maxT));
    maxSp = 1:length(maxT);
    d2 = [];
    Tmerged2 = Tmerged1;
%     if nargout > 1
%         [groupsSp maxSp d2] = ana.mergeSpikeTrains(gdf, Tmerged1, 'binSize', P.binSize,...
%             'minPercentOfEqualSPikes', P.minPercentOfEqualSPikes, ...
%             'nMaxSpikesPerSpikeTrain', P.nMaxSpikesPerSpikeTrain);
%     else
%         [groupsSp maxSp] = ana.mergeSpikeTrains(gdf, Tmerged1, 'binSize', P.binSize,...
%             'minPercentOfEqualSPikes', P.minPercentOfEqualSPikes, ...
%             'nMaxSpikesPerSpikeTrain', P.nMaxSpikesPerSpikeTrain);
%     end  
%     Tmerged2 = Tmerged1(:,:,maxSp);
    
    %% Add both indices
    tidx = maxT(maxSp);

    
    %%
    if nargout > 1
        D.templates.D = d1;
        D.templates.groups = groupsT;
        D.templates.maxT = maxT;
        D.templates.T = Tmerged1;
        D.spiketrains.D = d2;
        D.spiketrains.groups = groupsSp;
        D.spiketrains.maxT = maxSp;
        D.spiketrains.T = Tmerged2;
    end