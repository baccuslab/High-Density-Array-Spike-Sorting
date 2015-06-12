
function [C degrees_of_freedom] = diagonalSubspaceLoading(C, condNumber, debugflag)
    % inputs:
    %    C - Covariance matrix
    %    condNumber - minimal condition number, C should have
    % 
    % algorithm:
    %    the condition number of C is the biggest singular value (s1) of C
    %    devided by the smallest (sN): cond(C) = s1/sN
    %    to set the condition number of C to c we need to increase the
    %    smallest singular value so that s1/(sN+x) = c
    %    =>    s1 = c*sN + c*x     =>    x = (s1 - c*sN)/c = s1/c - sN 
    %    add so long the outer products of the eigenvectors of C to C
    %    whose eigenvalues are smaller than (sN+x).

    %    since C is assumed to be a covariance matrix, we can use the 
    %    eigenvalues instead of the singular values (they are the same)
    [V, D] = eig(C);
    D = diag(D);
    assert(~any(imag(D)), 'C has complex eigenvalues!');
    assert(~any(D<0), 'C has negative eigenvalues!');
    s1 = max(D);
    c  = condNumber;
    minEV = s1/c;
    if nargin>2
        mysort.plot.figure;
        plot(abs(D), '.-b'); hold on
    end
    degrees_of_freedom = size(C,1);
    count = 1;
    while any(D < minEV*.999)
        idx = find(D < minEV,1);
        degrees_of_freedom = degrees_of_freedom -1;
        x = minEV - D(idx);
        % THIS STEP IS NECESSARY! Do not calculate C = C + x*V*V'
        % directly!
        C_ = V(:,idx)*V(:,idx)';
        % OTHERWISE, C will not be symmetric anymore due to numerical 
        % instabilities and can get complex eigenvalues!
        C = C + x*C_;
        [V, D] = eig(C);
        D = diag(D);        
        if nargin>2
            if any(imag(D))
                fprintf('Warning in DSL: C has complex eigenvalues i=%d!\n',i);
            end
            plot(1:length(D), abs(D), '.-', 'color', mysort.plot.vectorColor(count));
            count = count+1;
        end
    end
end    