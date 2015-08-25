function [xticks Y] = templateSpikeTrainData(T, gdf, cutLeft)
    Tf = size(T,1);
    nC = size(T,2);
    nT = size(T,3);
    
    if nargin >4
        gdf = gdf(gdf(:,2)>=t1 & gdf(:,2)<=t2,:);
    end
    if isempty(gdf)
        xticks = [];
        Y = [];
        return
    end
    
    % add amplitudes if necessary
    if size(gdf,2) == 2
        gdf = [gdf ones(size(gdf,1),1)];
    end
    
    epochs = [gdf(:,2)-cutLeft gdf(:,2)-cutLeft + Tf-1];
    epochsM = mysort.epoch.merge(epochs);
    epochsMLength = mysort.epoch.length(epochsM);
    
    dataL = sum(epochsMLength) + size(epochsM,1);
    epochsMstartInData = cumsum(epochsMLength+1)';
    epochsMstartInData = [0 epochsMstartInData(1:end-1)]+1;
    
    xticks = zeros(1, dataL);
    
    for i=1:size(epochsM,1)
        s1 = epochsMstartInData(i);
        s2 = s1+epochsMLength(i)-1;
        xticks(s1:s2+1) = [epochsM(i,1):epochsM(i,2) nan];
    end
    
    Y = zeros(nC, dataL);
    
    myEpochIdx = 1;
    for i=1:size(gdf,1)
        if gdf(i,2) > epochsM(myEpochIdx,2)
            myEpochIdx = myEpochIdx +1;
        end
        
        offsetInEpoch = gdf(i,2)-cutLeft - epochsM(myEpochIdx,1);
        
        s1 = epochsMstartInData(myEpochIdx)+offsetInEpoch;
        s2 = s1+Tf-1;
        if s1>0 && s2 <= size(Y,2)
            Y(:, s1:s2) = Y(:, s1:s2)+gdf(i,3)*[squeeze(T(:,:,gdf(i,1)))'];
        end
%             warning()
    end
    
    
    
    

