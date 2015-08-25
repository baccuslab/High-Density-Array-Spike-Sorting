%%
allPeaks = sortrows(cell2mat(peaksOvp'));
keepIdx = zeros(size(allPeaks,1),1);
% find spikes
count = 1;
while count <=size(allPeaks,1)
    k = count+1;
    while (k<=size(allPeaks,1)) && ( (allPeaks(k,1) - allPeaks(count,1)) < 10)
        k=k+1;
    end
    [mx mxidx] = max(allPeaks(count:k-1,2));
    keepIdx(count+mxidx-1) = 1;
    count=k;
end
allPeaks = allPeaks(keepIdx==1,:);