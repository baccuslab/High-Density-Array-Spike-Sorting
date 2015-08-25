function Cce = Cte2Cce(Cte, nC)
% Compute the channel embedded (ce) representation of the noise block covariance
% matrix C with nC channels from the time embedded (te) representation.
% The time embedding is the one chosen e.g. in Pouzat2002 or Franke2010.
% The channel embedding is chosen more in the signal processing literature,
% e.g. in Luetkepohl2005 and is more numerically robust.

% Tf = size(Cte,1)/nC;

Tf = size(Cte, 1)/nC;
maxLag = Tf-1;
Cce = zeros(Tf*nC, Tf*nC);

for lag=0:maxLag
    ccollag = mysort.noise.Cte2CcolLag(Cte, nC, lag);
    if lag == 0
        for lag2 = 0:maxLag
            idx = lag2*nC+1 : (lag2+1)*nC;
            Cce(idx,idx) = ccollag;
        end
    else
        for lag2 = 0:maxLag-lag
            idx1 = (lag2+lag)*nC+1 : (lag2+lag+1)*nC;
            idx2 = lag2*nC+1 : (lag2+1)*nC;
            Cce(idx1,idx2) = ccollag;
            Cce(idx2,idx1) = ccollag';
        end
    end  
end
