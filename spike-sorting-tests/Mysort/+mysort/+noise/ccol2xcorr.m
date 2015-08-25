function xc = ccol2xcorr(Ccol)
    nC = size(Ccol,2);
    Tf = size(Ccol,1)/nC;
    maxLag = Tf-1;
    
    xc = zeros(2*maxLag+1, nC*nC);
    
    ccount = 1;
    for c1 = 1:nC
        for c2 = 1:nC
            for lag = -maxLag:maxLag
                if lag < 0 
                    xc(lag+maxLag+1, ccount) = Ccol(-lag*nC + c2, c1);
                else
                    xc(lag+maxLag+1, ccount) = Ccol( lag*nC + c1, c2);
                end
            end
            ccount = ccount+1;
        end
    end