function matches = findNearSpikes(st1, st2, maxDist)
    % Finds near spikes in st1 and st2 that are maximally maxDist apart.
    % Input:
    %    st1   - spike train 1
    %    st2   - spike train 2
    %  maxDist - maximal allowed distance between to spikes so that they
    %            are marked as being near
    % Output:
    %  matches - matrix with 3 columns. Every row corresponds to a spike
    %            pair between a spike in st1 and st2. First column is the
    %            index i1 into st1, second the index i2 into st2 and the
    %            third is the time distance of the two spikes as 
    %               dt = st2(i2) - st1(i1);
    %            every spike can appear with multiple partners in that list
    % Algorithm:
    %    The problem is symmetric in st1,st2.
    %    We arbitrarily choose st1 as the reference spike train.
    %    One pointer pRef points on spikes in st1, another pointer pFail 
    %    points on the last spike in st2 that could not be connected to a
    %    spike in st1.
    %    This is where we have to start looking for partners of our 
    %    reference spike in st1. The third pointer pMove is moved from here
    %    onwards. 
    %    As long as spikes cannot be matched the pMove and pFail are moved
    %    forward. If spikes are matched only pMove is moved until the first
    %    non matching spike is found. Then pRef is moved.
    if isempty(st1) || isempty(st2)
        return
    end
    n1 = length(st1);
    n2 = length(st2);
    
    % initialise with "guess" of the size
    matches = zeros(min(n1, n2)*100,3);
    
    pRef = 1;
    pFail = 0;
    pMove = [];
    nM = 1;
    while pRef <= n1
        pMove = pFail+1;
        while pMove <= n2
            dt = st2(pMove) - st1(pRef);
            if dt < -maxDist % check if move points before ref
                % fail before
                pFail = pFail +1;
                pMove = pMove +1;
            else
                % match of fail behind
                break
            end
        end
        
        % This seems to be the fastest way. Try catch is too slow
        % as well as precomputing the maximal move
        while pMove <= n2
            dt = st2(pMove) - st1(pRef);
            if abs(dt) <= maxDist
                % match
                matches(nM, :) = [pRef pMove dt];
                nM = nM+1;                
            else
                % fail behind
                break
            end
            pMove = pMove+1;
        end
        pRef = pRef+1;
    end
    matches = matches(1:nM-1,:);