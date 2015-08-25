function CCol = Cte2Ccol(Cte, nC)
Tf = size(Cte, 1)/nC;
maxLag = Tf-1;
CCol = zeros(maxLag*nC, nC);

for lag=0:maxLag
    idx = lag*nC+1 : (lag+1)*nC;
    CCol(idx,:) = mysort.noise.Cte2CcolLag(Cte, nC, lag);
end

