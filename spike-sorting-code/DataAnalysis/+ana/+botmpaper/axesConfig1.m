function axesConfig1(ax, nC)
    set(ax, 'box', 'off');
    set(ax, 'xticklabel', []);
    set(ax, 'yticklabel', []);
    set(ax, 'xcolor', 'w');
    set(ax, 'ycolor', 'w');
    if nC == 1
        axis(ax, 'tight');
    else
        set(ax, 'ylim', [0 nC*2+2])
    end
    