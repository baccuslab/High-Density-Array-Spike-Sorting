function T = QuirogaMatchingMoreTemplatesHelper(T, ti)
    % Function that builds a new template set from the template set in T
    [nT, Tf] = size(T);
    
    for i=1:ti
        n = randperm(2,1);
        t = zeros(1, Tf);
        r = randperm(3,n);
        tt = T(r,:);
        for k=1:n
            if rand>.8
                alpha = .25*randn + 2;
            else
                alpha = max(.1, .1*randn + .5);
            end
            t = t + alpha*tt(k,:);
        end
        T = [T; t];
    end