function C = xcov2timeC(xcovs, maxLag)
    % Computes the time embedded (Pouzat 2002) noise covariance matrix
    % from a set of x-covariance function in xcovs with maximal lag maxlag.
    
    nC = size(xcovs,1);
    L = length(xcovs{1,1});
    maxLagCalculated = (L-1)/2;
    if nargin == 1
        maxLag = maxLagCalculated;
    end
    Tf = maxLag+1;
    C = zeros(nC*(maxLag+1), nC*(maxLag+1));
    for row=1:nC
        r1 = (row-1)*Tf +1;
        r2  = row*Tf;        
        xc = xcovs{row,row};
        C(r1:r2, r1:r2) = toeplitz(xc(Tf:2*Tf-1), xc(Tf:-1:1) );
        for col=row+1:nC
            c1 = (col-1)*Tf +1;
            c2  = col*Tf;
            xc = xcovs{row,col};
            C(r1:r2, c1:c2) = toeplitz(xc(Tf:2*Tf-1), xc(Tf:-1:1) );
            C(c1:c2, r1:r2) = toeplitz(xc(Tf:2*Tf-1), xc(Tf:-1:1) )';
        end
    end