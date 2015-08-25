function [S a b ci groups] = matrixSortRowsColumns(M, method, parameter)
    % sorts the rows and columns of matrix M with a specific method
    % e.g. is ArgMax is selected, M will be resorted to S so that
    % S(1,1) = max(M) and S(2,2) = nextMax(M)
    % Input: 
    %   M         -  
    %   method    -
    %   parameter -
    %
    % Output:
    %   S         - resorted matrix
    %   a         -
    %   b         -
    %   ci        - only for method = 'components'
    %   groups    - only for method = 'components'
    if nargin < 2
        method = 'ArgMax';
    end
    S = M;
    ci = []; groups = {};
    if strcmp(method, 'ArgMax')
        n = min(size(M));
        a = zeros(n, 1);
        b = zeros(n, 1);
        
        for i=1:min(n)
            [a(i) b(i)] = mysort.util.matrixArgMax(M);
            M(a(i),:) = -inf;
            M(:,b(i)) = -inf;
        end
    elseif strcmp(method, 'components')
        assert(nargin == 3, 'If components is used, we need a sparse mask for M')
        A = parameter;
        ci = components(A);  
        uGroups = unique(ci);
        groups = cell(length(uGroups),1);
        for i=1:length(uGroups)
            groups{i} = find(ci==uGroups(i))';
        end
        a = cell2mat(groups);
        b = a;
    else
        error('unknown method: %s', method)
    end
    % add rows that were not associated (only if size(M,1) > size(M,2))
    a = [a; setdiff(1:size(M,1), a)'];
    % add cols that were not associated (only if size(M,1) < size(M,2))
    b = [b; setdiff(1:size(M,2), b)'];
    S = S(a,b);    
