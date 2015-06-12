function CCol = xcorr2CCol(xc)
%     maxlag = (size(xc,1)-1)/2;
%     Tf = maxlag+1;
    nC = sqrt(size(xc,2));
    Cte = mysort.noise.xcorr2Cte(xc);
    CCol = mysort.noise.Cte2Ccol(Cte, nC);