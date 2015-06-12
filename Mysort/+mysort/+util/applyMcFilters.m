
function Y = applyMcFilters(X, F)        
    nF = size(F,1);
    nC = size(X,1);
    Tf = size(F,2)/nC;
    Y = zeros(nF, size(X,2));
    for f=1:nF
        myF = F(f,:);
        Y(f,:) = mysort.util.mcfilt(X, mysort.wf.v2m(myF, nC));
        Y2(f,:) = mysort.util.mcfilt(X, mysort.wf.v2m(myF, nC), 'valid');
        Y3(f,:) = mysort.util.mcfilt(X, mysort.wf.v2m(myF, nC), 'full');
%         Y(f,1:floor(Tf/2)) = zeros(1,floor(Tf/2));
    end        
end