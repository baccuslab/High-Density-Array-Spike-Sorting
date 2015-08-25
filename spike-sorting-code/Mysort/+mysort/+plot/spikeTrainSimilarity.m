function P = spikeTrainSimilarity(R, varargin)
    P.fh = [];
    P.ah = [];
    P.sort = 1;
    P.normalize = 'gt'; % normalize the similarity (which is just number of common spikes) by the length of the "gt" spike trains
                        % alternative 'none' and 'sort' and both
    P.log = 0;          % perform a log transform of final plot
    P = mysort.util.parseInputs(P, varargin, 'error');
    
    if isempty(P.ah)
        if isempty(P.fh)
            P.fh = mysort.plot.figure('w', 800, 'h', 800);
        end
        P.ah = mysort.plot.axes();
    end
    if strcmp(P.normalize, 'both')
        if length(P.ah) ~= 2
            P.ah(1) = subplot(2,1,1);
            set(P.ah(1), 'fontSize', 14);
            P.ah(2) = subplot(2,1,2);
            set(P.ah(2), 'fontSize', 14);
        end
        P1 = P;
        P1.normalize = 'gt';
        P1.ah = P.ah(1);
        P1 = mysort.plot.spikeTrainSimilarity(R, P1);
        P2 = P;
        P2.normalize = 'sort';
        P2.ah = P.ah(2);
        P2.sort = P1.sort;
        mysort.plot.spikeTrainSimilarity(R, P2);
        linkaxes(P.ah, 'xy');
        return
    end
    S = R.similaritySt1St2;
    lastCol = R.nFN(:);
    lastRow = R.nFP2(:)';
    ylab = '# Spikes';
    if ~strcmp(P.normalize, 'none')
        % normalize FNs in gt with gt
        lastCol = lastCol./R.nSP1(:);
        % normalize FPs in sort with sort
        lastRow = lastRow./R.nSP2(:)';
        if strcmp(P.normalize, 'gt')
            % normalize S with gt
            S = S./repmat(R.nSP1(:), 1, size(S,2));
            ylab = '% Spikes Gt';
        elseif strcmp(P.normalize, 'sort')
            % normalize S with sort
            S = S./repmat(R.nSP2(:)', size(S,1), 1);
            ylab = '% Spikes Sort';
        else
            error('unknown normalization (either gt, sort, or none)');
        end
    end
    
    if ~isempty(P.sort) && (iscell(P.sort) || P.sort)
        if iscell(P.sort)
            a = P.sort{1};
            b = P.sort{2};
            S = S(a,b);
        else
            [S a b] = mysort.util.matrixSortRowsColumns(S);
        end
        lastCol = lastCol(a);
        lastRow = lastRow(b);
        P.sort = {a, b};
    end
    S = [S lastCol;
        lastRow 0];
    
    
    if P.log
        S = log(S);
        ylab = ['log(' ylab ')'];
    end
    
    imagesc(S, 'parent', P.ah);
    xlabel(P.ah, 'Sorted Neurons (last col FN)')
    ylabel(P.ah, 'Ground Truth (last row FP)');
    P.ch = colorbar('peer', P.ah);
%     set(get(P.ch, 'ylabel'), 'string', ylab);
    title(P.ah, ylab);
    