function st = mgdf2spikeTrains(mgdf)
    units = unique(mgdf(:,1));
    if size(mgdf,2) == 2
        mgdf(:,3) = 1;
    end
    
    trials = unique(mgdf(:,3));
    
    st = cell(length(trials), length(units));
    for ui=1:length(units)
        % Do this in two lines, not in one. Faster this way. (why?!)
        idx = mgdf(:,1)==units(ui);
        mgdf_ = mgdf(idx,:);
        for ti=1:length(trials)
            st{ti, ui} = mgdf_(mgdf_(:,3) == trials(ti), 2);
        end
    end
end