function [xc A S] = tXCorr(T)
    [Tf nC nT] = size(T);
    
    count = 1;
    xc = [];
    ij = [];
    
    for i=1:nT
        ti = squeeze(T(:,:,i));
        for j=i+1:nT
            tj = squeeze(T(:,:,j));
            xc(:,count) = xcov(ti(:,1), tj(:,1), 'biased');
            for c=2:nC
                xc(:,count) = xc(:,count) + xcov(ti(:,c), tj(:,c), 'biased');
            end
            xc(:,count) = xc(:,count)/(nC*std(ti(:))*std(tj(:)));
            ij(count,:) = [i j];
            count = count+1;
        end
    end
    
    if nargout > 1
        A = eye(nT);
        S = zeros(nT);
        [maxis shifts] = max(xc, [], 1);
        shifts = shifts - Tf;
        count = 1;
        for i=1:nT
            for j=i+1:nT
                A(i,j) = maxis(count);
                S(i,j) = shifts(count);
                A(j,i) = maxis(count);
                S(j,i) = -shifts(count);
                count = count +1;
            end
        end
    end
        