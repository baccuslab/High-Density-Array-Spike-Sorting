function CCol = xcov2ccol(xcovs, maxLag)
    % Computes the column representation of a noise covariance matrix from
    % x-covariance function in xcovs with maximal lag maxlag.
    
    nC = size(xcovs,1);
    L = length(xcovs{1,1});
    maxLagCalculated = (L-1)/2;
    if nargin == 1
        maxLag = maxLagCalculated;
    end
    
    CCol = zeros(nC*(maxLag+1), nC);
    tau0 = maxLagCalculated+1;
    for i=1:nC
        for j=i:nC
            if i==j
                % same channel, autocov is symmetric
                for tau = 0:min(maxLag, maxLagCalculated)
                    CCol(tau*nC + i,i) = xcovs{i,i}(tau0+tau);
                end
            elseif ~isempty(xcovs{i,j})
                % zero lag exists only once
                CCol(i, j) = xcovs{i,j}(tau0);
                CCol(j, i) = xcovs{i,j}(tau0);
                % all other lags exist + and -
                for tau = 1:min(maxLag, maxLagCalculated)
                    CCol(nC*tau + i, j) = xcovs{i,j}(tau0+tau);
                    CCol(nC*tau + j, i) = xcovs{i,j}(tau0-tau);
                end
            end
        end
    end