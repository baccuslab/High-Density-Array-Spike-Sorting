function Cte = xcorr2Cte(xc)
% This function computes the time embedded block covariance matrix from the
% xcorr functions in xc. xc is the matrix e.g. returned by the matlab xcorr
% function. Cte is the matrix used e.g. in Pouzat2002 or Franke2010.
maxlag = (size(xc,1)-1)/2;
Tf = maxlag+1;
nC = sqrt(size(xc,2));

Cte = zeros(Tf*nC, Tf*nC);

xccol = 1;
for i=1:nC
    iidx = (i-1)*Tf +1: i*Tf;
    for j=1:nC
        jidx = (j-1)*Tf +1: j*Tf;
        Cte(iidx, jidx) = toeplitz(xc(Tf:end, xccol), flipud(xc(1:Tf, xccol)));
        xccol = xccol +1;
    end
end