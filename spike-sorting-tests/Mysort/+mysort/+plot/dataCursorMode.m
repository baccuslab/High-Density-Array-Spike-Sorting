
function dataCursorMode(figHandle)
    persistent old_pos;
    persistent old_cval;
    if nargin == 0
        figHandle = gcf;
    end
    set(datacursormode(figHandle),'UpdateFcn',{@mydataCursor, '%10.1f'});

    % ---------------------------------------------------------------------
    function txt = mydataCursor(empt, event_obj, str)
        %global dtip;
        pos = get(event_obj,'Position');
        % check if this is an imagesc plot
        target = get(event_obj, 'Target');
        try
            cdata = get(target, 'CData');
            xdata = get(target, 'XData');
            ydata = get(target, 'YData');
            xidx = find(abs(xdata-pos(1))<eps);
            yidx = find(abs(ydata-pos(2))<eps);
            cval = cdata(yidx, xidx);
        catch
            cval = [];
        end
        
        px = sprintf(str,pos(1));
        py = sprintf(str,pos(2));
        txt = {['X: ', px],['Y: ', py]};
        if ~isempty(cval)
            txt = [txt {['C: ', sprintf(str, cval)]}];
        end
        dtip = struct();
        dtip.x = pos(1);
        dtip.y = pos(2);
        assignin('base','dtip',struct());
        fprintf('x: %10.2f y: %10.2f', dtip.x, dtip.y);
        if exist('old_pos','var')
            try
                dx = sprintf(str,pos(1)-old_pos(1));
                dy = sprintf(str,pos(2)-old_pos(2));
                txt = {txt{:}, ['dX: ', dx], ['dY: ', dy]};
                dtip.dx = pos(1)-old_pos(1);
                dtip.dy = pos(2)-old_pos(2);
                fprintf(' dx: %10.2f dy: %10.2f\n', dtip.dx, dtip.dy);
            catch
                % old_pos is not properly initialized. Ignore this.
                fprintf('\n');
            end
        else
            fprintf('\n');
        end
        old_pos = pos;
        old_cval = cval;
        dtip.old_pos = pos;
        assignin('base','dtip',dtip);
    end
end