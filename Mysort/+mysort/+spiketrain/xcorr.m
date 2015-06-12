
function [xc, E_xc, var_xc, conf_95_xc, bin_centers, edges, xc_unnormalized] = xcorr(st1, st2, varargin)
    % normalization in Dayan&Abbot is not good
    % Use the work from Brillinger 1976, Brillinger et al. 1976 and
    % Brillinger 1992. Makes much more sense.
    %
    % However, this here does only work with stationary processes. 
    %
    %   T            : time length
    %   Na(T), Nb(T) : number of spikes in whole spike train for neuron a
    %                  and neuron b
    %   nua, nub     : firing rates of neuron a  = Na(T)/T
    %   P(Nb(t) = n) : probability of n spikes in bin of length t for
    %                  neuron b   = (nub*t)^n /n!  * exp(-nub*t)
    %                  (follows from poisson process)
    %  
    P.maxLag  = 100;   % in samples
    P.binSize = 32;    % in samples
    P.T       = [];    % length of trial/spiketrain in sample
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
        E_xc    = zeros(nTrials, 1);
        var_xc  = zeros(nTrials, 1);
        conf_95_xc  = zeros(nTrials, 1);
        xc_unnormalized = xc;
        if length(P.T) == 1
            P.T = repmat(P.T, nTrials, 1);
        end
        not_null_trials = 0;
        for trial = 1:nTrials
            T = [];
            if ~isempty(P.T)
                T = P.T(trial);
            end
            if ~isempty(st1{trial}) && ~isempty(st2{trial})
                not_null_trials = not_null_trials+1;
                [xc(trial,:) E_xc(trial, :) var_xc(trial,:) conf_95_xc(trial,:) xc_unnormalized(trial,:)] = single_xcorr(st1{trial}, st2{trial}, T); 
            end
        end
        if not_null_trials>0
            xc              = sum(xc,              1)/not_null_trials;
            E_xc            = sum(E_xc,            1)/not_null_trials;
            var_xc          = sum(var_xc,          1)/not_null_trials;
            conf_95_xc      = sum(conf_95_xc,      1)/not_null_trials;
            xc_unnormalized = sum(xc_unnormalized, 1)/not_null_trials;
        else
            xc              = [];
            E_xc            = [];
            var_xc          = 1;
            conf_95_xc      = 2;
            xc_unnormalized = [];            
        end
    else
        [xc E_xc var_xc conf_95_xc xc_unnormalized] = single_xcorr(st1, st2, P.T);   
    end
    
    %----------------------------------------------------------------------
    function [xc E_xc var_xc conf_95_xc xc_unnorm] = single_xcorr(ST1, ST2, T)
        xc = zeros(1,nBins); xc_unnorm = xc;
        if isempty(ST1) || isempty(ST2)
            return
        end        
        if isempty(T)
            T = max(ST1(end), ST2(end)) - min(ST1(1), ST2(1));
            if T == 0
                T = max(max(ST1), max(ST2));
            end
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

        B  = P.binSize;
        m_A       = Na/T;
        m_AB      = xc/(Nb*B);
        m_B       = Nb/T;
        xc_unnorm = xc;
        xc        = m_AB;
        E_xc      = m_A;           % This is Brillingers m_A = N_A(T)/T
        var_xc    = Nb/(T*Na*B);
        conf_95_xc= 1/ (B*Nb);
    end
end