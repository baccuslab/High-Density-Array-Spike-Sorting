

function E = calculateEnergy(X, invC)
    % calculates (mahalanobis) energy of signal
    nC = size(X,1);
    Tf = size(invC,1)/nC;
    nS = size(X,2);
    E = zeros(1, nS);
    Tf2 = floor(Tf/2);
    for ii=1:nS-Tf
        Xt = X(:,ii:ii+Tf-1)';
        x = Xt(:);
        E(1,ii+Tf2) = x'*invC*x;
    end        
end  
