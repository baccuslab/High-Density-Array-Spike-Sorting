function plotDirectionSelectivity(ax, T, t, C, GC)
return
    e = load('U:\frameno_epochs_for_6March_3.mat');
    e = e.idx_where_to_cut;
    
    dirs = [270,90,0,180,225,45,315,135];
    reorder = [ 3 6 2 8 4 5 1 7 3];
    dirsre = dirs(reorder);
    t = T.templates(t);
    idx = t.getIdx();
    times = T.wfs.eventTimes(idx);
%     figure;
%     plot(times, zeros(length(times), 1), '.');
%     hold on
%     plot(e(:,1), zeros(size(e,1), 1), 'gx')
%     plot(e(:,2), zeros(size(e,1), 1), 'rx')
    
    for i=1:size(e,1)
        eidx = find(times >= e(i,1) & e(i,2) >= times);
        c(i) = length(eidx);
    end
%     meandir = sum(dirs.*c/sum(c));
    cla(ax);
    h = polar(ax, deg2rad(dirsre), c(reorder)/sum(c));
    set(h, 'color', 'b', 'linewidth', 2); 
%     set(ax,'nextplot', 'add');
%     h = polar(ax, [0 deg2rad(meandir)], [0 max(c)]);
%     set(h, 'color', 'k', 'linewidth', 3); 
    xlabel(ax,'orientation')
    ylabel(ax,'spike count')
    
    
    
    