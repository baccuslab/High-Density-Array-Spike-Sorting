function mgdf = mergeGdfs(gdfList)
    nGdfs = length(gdfList);
    nSpikesPerGdf = zeros(1, nGdfs);
    for i=1:length(gdfList)
        nSpikesPerGdf(i) = size(gdfList{i},1);
    end
    mgdf = zeros(sum(nSpikesPerGdf), 3);
    for i=1:length(gdfList)
        s1 = sum(nSpikesPerGdf(1:i-1))+1;
        s2 = sum(nSpikesPerGdf(1:i));
        
        mgdf(s1:s2,:) = [repmat(i, nSpikesPerGdf(i), 1) gdfList{i}];
    end