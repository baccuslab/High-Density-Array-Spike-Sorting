
function [xc bin_centers not_null_trials] = xcorr_raw(st1, st2, varargin)
    % raw xcorr in spike counts 
    P.maxLag  = 100;   % in samples
    P.binSize = 32;    % in samples
    P = mysort.util.parseInputs(P, 'xcorr', varargin);

    % center edges of histogram so that the middle bin has its center on
    % zero
    nBins_one_direction = floor(P.maxLag/P.binSize);
    nBins = 2*nBins_one_direction + 1;
    
    maxDist = nBins_one_direction * P.binSize + P.binSize/2;
    edges = linspace(-maxDist, maxDist, nBins+1);
    if nargout > 1
        bin_centers = edges(2:end)-P.binSize/2;
    end

    if iscell(st1)
        nTrials = length(st1);
        assert(nTrials == length(st2), 'number of trials in st1 and st2 must be identical!');
        xc      = zeros(nTrials, nBins);
        not_null_trials = 0;
        for trial = 1:nTrials
            if ~isempty(st1{trial}) && ~isempty(st2{trial})
                xc(trial,:) = single_xcorr(st1{trial}, st2{trial}); 
            end
        end
    else
        xc = single_xcorr(st1, st2);   
    end
    
    %----------------------------------------------------------------------
    function xc = single_xcorr(ST1, ST2)
        xc = zeros(1,nBins); 
        if isempty(ST1) || isempty(ST2)
            return
        end        
        Na = length(ST1);
        Nb = length(ST2);
        activeSt1     = 1;
        firstValidSt2 = 1;
        currentSt2    = 1;
        while activeSt1 <= Na
            currentDist = ST2(currentSt2) - ST1(activeSt1);
            firstValidFound = false;
            while currentDist < maxDist % had to change <= to < (see below)
                if currentDist >= -maxDist
                    if ~firstValidFound
                        firstValidFound = true;
                        firstValidSt2 = currentSt2;
                    end
                    % This works, but is slow and produces an empty idx for
                    % the case currentDist == maxDist
                    % binidx1 = find(edges>currentDist, 1) -1;

                    % This should be faster but made it necessary to ignore
                    % the case currentDist == maxDist
                    binidx = floor(((maxDist+currentDist)/P.binSize))+1;
                    % assert(binidx1 == binidx, 'binidx diff');
                    xc(binidx) = xc(binidx)+1;
                end

                currentSt2 = currentSt2 +1;
                if currentSt2 > Nb
                    break
                end
                currentDist = ST2(currentSt2) - ST1(activeSt1);
            end
            currentSt2 = firstValidSt2;
            activeSt1 = activeSt1+1;
        end
    end
end