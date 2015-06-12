function [dt1, dt2] = computeDeltaTauBetweenTwo(st1, st2)
    % dump quick and dirty implementation
    dt1 = zeros(1,length(st1));
    dt2 = zeros(1,length(st2));
    for i=1:length(st1)
        dt = st2-st1(i);
        [~, idx] = min(abs(dt));
        dt1(i) = dt(idx); 
    end
    for i=1:length(st2)
        dt = st1-st2(i);
        [~, idx] = min(abs(dt));
        dt2(i) = dt(idx); 
    end
end