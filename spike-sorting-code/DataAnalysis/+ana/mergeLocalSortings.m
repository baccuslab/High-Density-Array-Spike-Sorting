function [G] = mergeLocalSortings(G, meanNoiseStd, varargin)
    % Merges the individual neurons from many local sortings of the same
    % piece of data. Here, a given real neuron might have been sorted in
    % multiple local sortings. Only the sorted neuron from the local sorting 
    % where the neuron had the best SNR should be kept, the others discarded
    % Input:
    %  G                   - array of struct containing the local sortings
    %  G(i).nSpikesPerElectrode - number of spikes for each template and
    %                        electrode that were averaged together to 
    %                        compute the template
    %  meanNoiseStd        - the average noise standard deviation
    % Output:
    %  groups   - groupings of templates. Each member in that group is
    %             thought to partly reflect the same neuron. Only one of
    %             those templates "survives" and is actually considered
    %             further
    %  maxT     - For each group in "groups" this is the index of the
    %             template in that group that is considered to be the best
    %             representative of that group
    %  D        - Debugging stuff
    %
    P.minSpikeTrainOvp = .50;       % minimal fraction of identical spikes
    P.overlapDist = 10;
    P.minCorrelationForMerger = .85; 
    P.minCorrelationForSureMerger = .95;
    P.percentageCutoff = .3;
    P.maxPercentageAmplitudeDistance = .6;
    P.sameHighestElectrodesPercent = 80;  % this is the percentage of valid electrodes on which both tempaltes have the highest energy
    P = mysort.util.parseInputs(P, varargin, 'error');
    
    % First step, walk through all templates and check if the maximal
    % amplitude is on an electrode that was inside the local sorting from
    % where the template came
    
    for g=1:length(G)
        if isempty(G(g).templates.wfs)
            nTg = 0; % size returns 1 for size(x,dim) with dim > ndims(x) which means size([], 3) == 1
            G(g).templates.duplicateTemplates = {};
            G(g).templates.masterTemplate = {};
            G(g).templates.maxInThisGroup = [];
        else
            nTg = size(G(g).templates.wfs,3);
            G(g).templates.maxInThisGroup = zeros(1,nTg);
            [mi ma mi_idx ma_idx] = mysort.wf.tMinMaxPerTemplate(G(g).templates.wfs);
            [globMin minChan] = min(mi, [], 2);

            for t=1:nTg
                G(g).templates.duplicateTemplates{t} = [];
                G(g).templates.masterTemplate{t} = [];            
                if ismember(minChan(t), G(g).sortedElectrodes)
                    G(g).templates.maxInThisGroup(t) = 1;
                end
            end            
        end        
    end
    
    % Now we know which template had his maximum actually in the group
    % where it was found. Let's assume that each neuron is also at least
    % found by the group of electrodes where its maximum was. Then, we can
    % safely ignore all templates from now on, which did not have that.
    %
    % However, the others might still be duplicates since a given electrode
    % can be in several groups and thus a neuron might be found in several
    % groups that actually had the electrode with the maximal amplitude.
    % So now we have to cross check all groups that share electrodes and
    % check if there are duplicates.  
    
    for g1=1:length(G)
        if isempty(G(g1).templates.wfs)
            nTg = 0;
        else
            nTg = size(G(g1).templates.wfs,3);
        end
        thisGroupElectrodes = G(g1).sortedElectrodes;
        if nTg<1
            continue
        end
        t1idx = find(G(g1).templates.maxInThisGroup);
        if isempty(t1idx)
            continue
        end
        for g2 = g1+1:length(G)
            if size(G(g2).templates.wfs,3)<1
                continue
            end
            t2idx = find(G(g2).templates.maxInThisGroup);
            if isempty(t2idx)
                continue
            end    
%             [C ia ib] = intersect(thisGroupElectrodes, G(g2).sortedElectrodes);
%             if isempty(C)
%                 continue
%             end
            % we found a pair with at least one electrode of overlap and
            % both groups have templates that had their maximal amplitude
            % in the respective group. Make sure those are not the same
            %
            % we do that by taking the biggest channels that have at
            % least a peak of 10% from the maximal amplitude and compute on
            % those the correlation between the maxima and minima of the            
            [mi1 ma1 mi_idx1 ma_idx1] = mysort.wf.tMinMaxPerTemplate(G(g1).templates.wfs(:,:,t1idx));
            [mi2 ma2 mi_idx2 ma_idx2] = mysort.wf.tMinMaxPerTemplate(G(g2).templates.wfs(:,:,t2idx));
            if 0
                figure;
                plot(mysort.wf.t2v(G(g1).templates.wfs(:,:,:))', 'b');
                hold on
                plot(mysort.wf.t2v(G(g2).templates.wfs(:,:,:))', 'g');
                figure;
                plot(mysort.wf.t2v(G(g1).templates.wfs(:,:,2))');
                hold on
                plot(mysort.wf.t2v(G(g2).templates.wfs(:,:,2))', 'g');
            end
            for t1 = 1:length(t1idx)
                M1 = abs([mi1(t1,:) ma1(t1,:)]);
                [sortedAmps1 sortedElectrodeIdices1] = sort(M1);
                absMax1 = max(M1);
                validChanIdx1 = M1>P.percentageCutoff*absMax1;
                validSortedIdx1 = sortedElectrodeIdices1(sortedAmps1>P.percentageCutoff*absMax1);
                for t2 = 1:length(t2idx)
                    bMerge = 0;
                    if 0
                        figure; plot(mysort.wf.t2v(G(g1).templates.wfs(:,:,:))');
                        figure;
                        plot(mysort.wf.t2v(G(g1).templates.wfs(:,:,t1))');
                        hold on
                        plot(mysort.wf.t2v(G(g2).templates.wfs(:,:,t2))', 'g');
                    end                    
                    M2 = abs([mi2(t2,:) ma2(t2,:)]);
                    absMax2 = max(M2);
                    validChanIdx2 = M2>P.percentageCutoff*absMax2;
                    validChanIdx = validChanIdx1 & validChanIdx2;
                    if (g1 == 8) & (g2 ==11) & (t1idx(t1)-1 == 5) & t2idx(t2)-1==2
                        disp('bla')
                    end
                    if ~any(validChanIdx)
                        continue
                    end
                    [sortedAmps2 sortedElectrodeIdices2] = sort(M2);
                    validSortedIdx2 = sortedElectrodeIdices2(sortedAmps2>P.percentageCutoff*absMax2);
                    highestValidElectrodeIntersect = length(intersect(validSortedIdx1, validSortedIdx2))/min(length(validSortedIdx1),length(validSortedIdx2));
                    if 100*highestValidElectrodeIntersect < P.sameHighestElectrodesPercent
                        continue
                    end
                    
                    MM1 = M1(validChanIdx);
                    MM2 = M2(validChanIdx);
                    c = MM1*MM2'/(norm(MM1)*norm(MM2));
                    if c > P.minCorrelationForMerger
                        if c > P.minCorrelationForSureMerger
                            bMerge = 1;
                            ovp = -1;
                            maxPercDist = -1;
                        else
                            % Templates are correlated among channels
                            [maxPercDist maxDistIdx] = max(abs(MM1-MM2));
                            maxPercDist = maxPercDist / max(abs(MM1(maxDistIdx)), abs(MM2(maxDistIdx)));
                            if maxPercDist < P.maxPercentageAmplitudeDistance
                                % Templates do not have a big difference
                                st1 = G(g1).gdf(G(g1).gdf(:,1) == t1idx(t1)-1,2);
                                st2 = G(g2).gdf(G(g2).gdf(:,1) == t2idx(t2)-1,2);
                                [O nO] = mysort.spiketrain.checkForOverlaps({st1, st2}, P.overlapDist);
                                ovp = nO(1)/min([length(st1) length(st2)]);
                                if c > P.minCorrelationForSureMerger || ovp > P.minSpikeTrainOvp
                                    % Spike trains are similar enough
                                    if 0
                                        mysort.plot.templates2D(G(g1).templates.wfs(:,:,t1idx(t1)), G(g1).templates.MES.electrodePositions, 100, 5, 'IDs', 1:length(t1idx))
                                        ah = gca; 
                                        hold on
                                        mysort.plot.templates2D(G(g2).templates.wfs(:,:,t2idx(t2)), G(g2).templates.MES.electrodePositions, 100, 5, 'IDs', length(t1idx)+(1:length(t2idx)), 'ah', ah)
                                    end
                                    bMerge = 1;
                                else
                                    % Rejection because of spike trains
                                    bMerge = 0;
                                    if 0
                                        mysort.plot.templates2D(G(g1).templates.wfs(:,:,t1idx(t1)), G(g1).templates.MES.electrodePositions, 100, 5, 'IDs', 1:length(t1idx))
                                        ah = gca; 
                                        hold on
                                        mysort.plot.templates2D(G(g2).templates.wfs(:,:,t2idx(t2)), G(g2).templates.MES.electrodePositions, 100, 5, 'IDs', length(t1idx)+(1:length(t2idx)), 'ah', ah)
                                    end                                
                                end
                            end
                        end
                    end
                    if bMerge
                        if sum(abs(M1)) > sum(abs(M2))
                            G(g1).templates.duplicateTemplates{t1idx(t1)} = [G(g1).templates.duplicateTemplates{t1idx(t1)}; [g2 t2idx(t2) c maxPercDist ovp]];
                            G(g2).templates.masterTemplate{t2idx(t2)} = [G(g2).templates.masterTemplate{t2idx(t2)}; [g1 t1idx(t1) c maxPercDist ovp]];
                        else
                            G(g2).templates.duplicateTemplates{t2idx(t2)} = [G(g2).templates.duplicateTemplates{t2idx(t2)}; [g1 t1idx(t1) c maxPercDist ovp]];
                            G(g1).templates.masterTemplate{t1idx(t1)} = [G(g1).templates.masterTemplate{t1idx(t1)}; [g2 t2idx(t2) c maxPercDist ovp]];
                        end
                    end
                end
            end
            if 0
                mysort.plot.templates2D(G(g1).templates.wfs(:,:,t1idx), G(g1).templates.MES.electrodePositions, 10, 5, 'IDs', 1:length(t1idx))
                ah = gca; 
                hold on
                mysort.plot.templates2D(G(g2).templates.wfs(:,:,t2idx), G(g2).templates.MES.electrodePositions, 10, 5, 'IDs', length(t1idx)+(1:length(t2idx)), 'ah', ah)
            end
        end
    end
        
    % Now, every template is marked as either
    % - not having its maximal amplitude in its local sorting
    % - having a duplicate template that is the same but smaller
    % - having a master template that is the same but bigger