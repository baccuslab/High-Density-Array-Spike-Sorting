function [msgdf nS] = gdf2multiSessionGdf(gdf, start_end_times, sessionidx, assert_flag)
    % Inputs:
    %    start_end_times    - can either be a matrix with two columns or a
    %                         vector
    %                         if vector:
    %                           is interpreted as consecutive periods where
    %                           the end of one period is the start of the
    %                           other. Periods connot be overlapping and
    %                           cover a consecutive piece of time.
    %                         if matrix:
    %                           first column is start time, second end time
    %                           for each period. Periods can be overlapping
    %                           or not completely covering

    % if the individual sesssion lengths are in L then
    % start_end_times = [0 cumsum(L)];
    if isempty(gdf)
        return
    end
    if nargin < 4
        assert_flag = 1;
    end
    
    if any(size(start_end_times) == 1)
        s1 = start_end_times(1:end-1);
        s2 = start_end_times(2:end);
        start_end_times = [s1(:) s2(:)];        
    end
    
    nS = size(start_end_times, 1);
    if isempty(start_end_times) || nS == 0
        error('star_end_times invalid');
    end
    if nargin == 2 || isempty(sessionidx)
        sessionidx = 1:nS;
    end
    msgdf = [gdf zeros(size(gdf,1),1)];
    for i=1:nS
        idx = gdf(:,2) > start_end_times(i,1) & gdf(:,2) <= start_end_times(i,2);
        msgdf(idx,3) = sessionidx(i);
        msgdf(idx,2) = msgdf(idx,2) - start_end_times(i,1);
        nS(i) = sum(idx);
    end
    assert(~any(msgdf(:,3)==0) || ~assert_flag, 'Some events were in no session!');