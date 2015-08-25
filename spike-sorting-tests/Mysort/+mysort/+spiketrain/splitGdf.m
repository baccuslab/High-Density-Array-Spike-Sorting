function gdfList = splitGdf(gdf, start_end_times, sessionidx)
    gdfList = {};
    if isempty(gdf)
        return
    end
    nS = length(start_end_times)-1;
    if isempty(start_end_times) || nS == 0
        gdfList{1} = gdf;
        return
    end
    
    if nargin < 3
        sessionidx = 1:nS;
    end
    for i=1:nS
        % check if this only one spike train
        if size(gdf,2) == 1
            tmpgdf = gdf(gdf > start_end_times(i) & gdf <= start_end_times(i+1),:);
            tmpgdf = tmpgdf - start_end_times(i);
            gdfList{sessionidx(i)} = tmpgdf;
        else
            tmpgdf = gdf(gdf(:,2) > start_end_times(i) & gdf(:,2) <= start_end_times(i+1),:);
            tmpgdf(:,2) = tmpgdf(:,2) - start_end_times(i);
            gdfList{sessionidx(i)} = tmpgdf;
        end
    end