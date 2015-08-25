
function ax = twoYAxis(x, y1, y2, xl, yl1, yl2, args1, args2)
    if ~exist('args1', 'var')
        args1 = {};
    end
    if ~exist('args2', 'var')
        args2 = {};
    end
    
    ax(1) = axes('yaxislocation', 'right'); box off
    pos = get(ax(1), 'position');
    plot(x, y2);
    ylabel(yl2);
    set(ax(1),'yaxislocation', 'right')
    
    
    ax(2) = axes('position', pos);
    plot(x, y1, args1{:});
    ylabel(yl1)
    xlabel(xl);

    ax(3) = axes('position', pos, 'yaxislocation', 'right'); box off
    plot(x, y2, args2{:});
    set(ax(3), 'visible', 'off');

   