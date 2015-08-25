function [gdf_resolved gdf_raw] = DiscriminantFunctions2Gdf(Y, Dmax, IDs, OvpIndex, refactory, maxTau, prior)
PRIOR = log(.99);
nT = size(Dmax, 2);
peaks  = {};
for i=1:nT
    [pks, locs] = findpeaks(Dmax(:,i), 'MINPEAKDISTANCE', refactory, 'MINPEAKHEIGHT', log(.99));
    % Check that 
    % a) a peak is also maximal over all other functions
    [mx mxidx] = max(Dmax(locs,:), [], 2);
    pks  = pks(mxidx == i);
    locs = locs(mxidx == i);
    ids  = IDs(locs, i);
    
    % b) that it is indeed the highest value in a short window around it
    idx = repmat(-refactory:refactory, length(locs), 1) + repmat(locs(:), 1, 2*refactory+1);
    idx = idx';
    idx(idx<1) = 1;
    idx(idx>size(Dmax,1)) = size(Dmax,1);
    
    X = reshape(Dmax(idx(:), i), 2*refactory+1, length(locs));
    [mx mxidx] = max(X,[], 1);
    keep = mxidx == refactory+1;
    pks  = pks(keep);
    locs = locs(keep);
    ids  = ids(keep);
    
    % c) if it is one of the boarder overlap functions, look for a peak in
    % that direction of an individual spike
    maxTauUsed = max(OvpIndex(:,4));
    myLeftBoarderOvpIDs = OvpIndex( abs(OvpIndex(:,4) + maxTauUsed)<1  &  OvpIndex(:,2)==i, 1);
    boarderPeaks = find(ismember(ids,myLeftBoarderOvpIDs));
    IDX = (repmat(-refactory:0, length(boarderPeaks), 1) + repmat(locs(boarderPeaks(:)), 1, refactory+1))';
    X = reshape(Y(IDX(:),i), refactory+1, length(boarderPeaks));
    [mx mxidx] = max(X,[], 1);
    keep = find( mx > PRIOR); %mxidx > 1 & mxidx < refactory &
    addpks = mx(keep)';
    addlocs = mxidx(keep)' + locs(boarderPeaks(keep))-refactory;
    addids = repmat(i, length(keep),1);
    pks(boarderPeaks) = []; 
    locs(boarderPeaks) = []; 
    ids(boarderPeaks) = []; 
    
    % SO FAR WE ONLY DELETED BOARDER PEAKS
    
    pks  = [pks; addpks];
    locs = [locs; addlocs];
    ids  = [ids; addids];
    
    % CHECK THE OTHER SIDE
    
    myRightBoarderOvpIDs = OvpIndex( abs(OvpIndex(:,4) - maxTauUsed)<1  &  OvpIndex(:,2)==i, 1);
    boarderPeaks = find(ismember(ids,myRightBoarderOvpIDs));
    IDX = (repmat(0:refactory, length(boarderPeaks), 1) + repmat(locs(boarderPeaks(:)), 1, refactory+1))';
    X = reshape(Y(IDX(:),i), refactory+1, length(boarderPeaks));
    [mx mxidx] = max(X,[], 1);
    keep = find( mx > PRIOR); %mxidx > 1 & mxidx < refactory &
    addpks = mx(keep)';
    addlocs = mxidx(keep)' + locs(boarderPeaks(keep));
    addids = repmat(i, length(keep),1);
    pks(boarderPeaks) = []; 
    locs(boarderPeaks) = []; 
    ids(boarderPeaks) = []; 
    
    pks  = [pks; addpks];
    locs = [locs; addlocs];
    ids  = [ids; addids];    
    
    peaks{i,1} = [locs pks ids];
end

allPeaks = sortrows(cell2mat(peaks));

% b) there is no bigger peak nearby from another function
keepIdx = zeros(size(allPeaks,1),1);
% find spikes
count = 1;
while count <=size(allPeaks,1)
    k = count+1;
    while (k<=size(allPeaks,1)) && ( (allPeaks(k,1) - allPeaks(count,1)) < maxTau)
        k=k+1;
    end
    [mx mxidx] = max(allPeaks(count:k-1,2));
    keepIdx(count+mxidx-1) = 1;
    count=k;
end
allPeaks = allPeaks(keepIdx==1,:);

    
% c) if is is an overlapfunction, there should be a peak from the other
% overlap function
keepIdx = ones(size(allPeaks,1),1);
% find spikes
k = 0;
while k <size(allPeaks,1)
    k = k+1;
    if keepIdx(k)==0
        continue
    end
    if allPeaks(k, 3) > nT
        % This is an overlapping spike
        myTime = allPeaks(k, 1);
        if abs(myTime-430962) < 25
            disp('d');
            allPeaks(abs(allPeaks(:,1)-myTime)<25 ,3)
        end
        myAmp = allPeaks(k, 2);
        myId = allPeaks(k, 3);
        myOvpIndex = OvpIndex(OvpIndex(:,1)==myId,:);
        myShift = myOvpIndex(4);
        mySingleId = myOvpIndex(2); 
        myPartner = myOvpIndex(3);
        idx = find(abs( (allPeaks(k, 1)-allPeaks(:, 1)) - myShift) < 6);
        idx = idx(allPeaks(idx,3)~=myId);
        if isempty(idx)
            % there is no other spike like me... keep me or throw me away?
            keepIdx(k) = 1;
        else 
            theirTime = allPeaks(idx, 1);
            theirAmp = allPeaks(idx, 2);
            theirId = allPeaks(idx, 3);     
            dt = myTime-theirTime;
            for j=1:length(theirId)
                if theirId(j) > nT
                    % this is also an overlap
                    theirOvpIndex = OvpIndex(OvpIndex(:,1)==theirId(j),:);
                    if theirAmp(j) < myAmp
                        keepIdx(idx(j)) = 0;
                    else
                        keepIdx(k) = 0;
                    end
                else
                    % this is not an overlap. check if it is larger than we are
                    if theirAmp(j) > myAmp
                        % yes, we are useless
                        keepIdx(k) = 0;
                    else
                        % no. the other one is useless
                        keepIdx(idx(j)) = 0;
                    end
                end
            end
        end
    end
end
allPeaks = allPeaks(keepIdx==1,:);

%

gdf_raw = allPeaks(:,[3 1]);

[gdf_resolved wasResolved] = mysort.spiketrain.resolveOvpIndexInGdf(gdf_raw, OvpIndex);
