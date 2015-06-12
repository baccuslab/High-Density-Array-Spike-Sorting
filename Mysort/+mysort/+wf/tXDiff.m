function [xc A S maxAbsDist maxAbsRelDist] = tXDiff(T)
    [Tf nC nT] = size(T);
    
    count = 1;
    xc = [];
    ij = [];
    
    for i=1:nT
        ti = squeeze(T(:,:,i));
        Eti = sum(ti(:).^2);
        for j=i+1:nT
            tj = squeeze(T(:,:,j));
            Etj = sum(tj(:).^2);
            xc(:,count) = xcorr(ti(:,1), tj(:,1), 'none');
            for c=2:nC
                xc(:,count) = xc(:,count) + xcorr(ti(:,c), tj(:,c), 'none');
            end
            xc(:,count) = Eti - 2*xc(:,count) + Etj;
            ij(count,:) = [i j];
            count = count+1;
        end
    end
    
    if nargout > 1
        T = mysort.wf.t2v(T);
        A = zeros(nT);
        S = zeros(nT);
        maxAbsDist = zeros(nT);
        maxAbsRelDist = zeros(nT);
        [mins shifts] = min(xc, [], 1);
        shifts = shifts - Tf;
        count = 1;
        for i=1:nT
            ti = T(i,:);
            for j=i+1:nT
                tj = T(j,:);
%                 figure; plot([ti; tj]')
                if abs(shifts(count)) <= 2
                    tj = mysort.util.shiftMCRows(tj, shifts(count), nC,1);
                    t = mysort.wf.vSubsel([tj; ti], nC, abs(shifts(count))+1:Tf-abs(shifts(count)));
                else
                    t = [tj; ti];
                end
%                 figure; plot(t')
                A(i,j) = mins(count);
                S(i,j) = shifts(count);
                A(j,i) = mins(count);
                S(j,i) = -shifts(count);
                DIST = t(1,:)-t(2,:);
                [maxAbsDist(i,j) maxDistIdx] = max(abs(DIST));
                maxAbsRelDist(i,j) = maxAbsDist(i,j)/max(abs(t(1,maxDistIdx)), abs(t(2,maxDistIdx)));
                maxAbsDist(j,i) = maxAbsDist(i,j);
                maxAbsRelDist(j,i) = maxAbsRelDist(i,j);
                count = count +1;
            end
        end
        
        
    end
        