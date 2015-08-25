function msgdf = gdf2multiSessionGdf(gdf, start_end_times, sessionidx, assert_flag)
    % if the individual sesssion lengths are in L then
    % start_end_times = [0 cumsum(L)];
    if isempty(gdf)
        return
    end
    if nargin < 4
        assert_flag = 1;
    end
    nS = length(start_end_times)-1;
    if isempty(start_end_times) || nS == 0
        error('star_end_times invalid');
    end
    if nargin == 2 || isempty(sessionidx)
        sessionidx = 1:nS;
    end
    msgdf = [gdf zeros(size(gdf,1),1)];
    for i=1:nS
        idx = gdf(:,2) > start_end_times(i) & gdf(:,2) <= start_end_times(i+1);
        msgdf(idx,3) = sessionidx(i);
        msgdf(idx,2) = msgdf(idx,2) - start_end_times(i);
    end
    assert(~any(msgdf(:,3)==0) || ~assert_flag, 'Some events were in no session!');