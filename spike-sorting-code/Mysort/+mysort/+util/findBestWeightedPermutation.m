function S = findBestWeightedPermutation(Q, W)
    % Find best permutation of Q so that similar (according to Q and
    % weighted by W) columns are next to each other (circular)

    % find for every dimension the permutation that it would favour.
    % Weight the pairwise values with the weight vectors from W
    nS = size(Q,1);
    P = zeros(nS, nS-1);
    V = P;
    for i=1:nS
        idx = setdiff(1:nS, i);
        WR = [W(i,idx)' idx'];
        WR = sortrows(WR);
        CR = [Q(i,idx)' idx']; 
        CR = sortrows(CR);
        V(i,:) = (WR(:,1).*CR(:,1))';
        P(i,:) = CR(:,2)';
    end
    % now V contains the values we get if we put a next to b. But if we
    % do that, we will also put b next to a. So we need to add this
    % value directly to the a-b value and vice versa:
    K = zeros(size(V));
    for a=1:nS
        for abidx = 1:nS-1
            b = P(a, abidx);
            baidx = find(P(b,:) == a);
            vab = V(a, abidx);
            vba = V(b, baidx);
            K(a, abidx) = vab+vba;
            K(b, baidx) = vab+vba;
        end
    end

    % walk along the list of preferred tupels and greedily accept 
    % them as long as S is not full 
    isBlocked = zeros(nS, 1);
    done = 0;
    J = K;
    % accept first tupel
    [a bidx] = mysort.util.matrixArgMax(J);
    J(a,bidx) = -inf;
    b = P(a, bidx);
    S = [a b zeros(1, nS-2)];
    pendingTuples = [];
    iter = 0;
    while ~done && iter < 100
        [a bidx] = mysort.util.matrixArgMax(J);
        J(a,bidx) = -inf;
        b = P(a, bidx);
        changedS = 0;
        sa = find(S==a,1);
        sb = find(S==b,1);          
        if isPossibleTupel(isBlocked, a, b, sa, sb)
            addTupel;
            if changedS==0
                addPendingTupel()
            end
        end
        
        if weAreDone(); return; end
        % if S was not changed pending tuples now fit
        
        if changedS
            checkPendingTupel()
            if weAreDone(); return; end
        end
        
        iter = iter+1;
    end
    % ---------------------------------------------------------------------
    function bo = isPossibleTupel(isBlocked, a, b, sa, sb)
        oneIsBlocked = isBlocked(a) || isBlocked(b);
        bothInList = ~isempty(sa) && ~isempty(sb);
        bo = ~oneIsBlocked && ~bothInList;
    end    
    % ---------------------------------------------------------------------
    function addTupel
        % check if a is already in S, then add b to its free side
        if ~isempty(sa)
            % a is already placed
            isBlocked(a) = 1;
            sl = mod(sa-1-1, nS)+1;
            sr = mod(sa-1+1, nS)+1;                    
            if S(sl)>0
                S(sr) = b;
            else
                S(sl) = b;
            end
            changedS = 1;
        elseif ~isempty(sb)
            % b is already placed
            isBlocked(b) = 1;
            sl = mod(sb-1-1, nS)+1;
            sr = mod(sb-1+1, nS)+1;                    
            if S(sl)>0
                S(sr) = a;
            else
                S(sl) = a;
            end          
            changedS = 1;
        end      
    end
    % ---------------------------------------------------------------------
    function checkPendingTupel()
        removeTupel = [];
        for ptup=1:size(pendingTuples)
            a = pendingTuples(ptup,1);
            b = pendingTuples(ptup,2);
            sa = find(S==a,1);
            sb = find(S==b,1);              
            if isPossibleTupel(isBlocked, a, b, sa, sb)
                addTupel;
            else
                removeTupel(end+1) = ptup;
            end
        end    
        pendingTuples(removeTupel,:) = [];
    end
    % ---------------------------------------------------------------------
    function bo = weAreDone()
        bo = ~any(S==0);
    end
    % ---------------------------------------------------------------------
    function addPendingTupel()
        % none is placed, add as pending tuple and try to add
        % in next iteration
        pendingTuples(end+1,:) = [a b];
    end
end