
function MAHA = calculateMahalanobisDistances(E, Y, M)       
    % Mahalanobis Distance
    nF = size(Y,1);
    nS = size(E,2);
    Tf = floor(size(M,1)/2);
    MAHA = zeros(nF+(nF*(nF-1)/2), nS);
    for f=1:nF
        MAHA(f,:) = E-2*Y(f,:) + M(Tf,f,f);
    end
    % Set negative mahalanobis values to zero.
    MAHA(MAHA<0) = 0;  
end