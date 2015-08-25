function F = toBins(st, edges)
    % this funny (x(:)')' is necessary if there is an empty spike train
    F = cellfun(@(x) st2f(x(:)')', st, 'uniformOutput', false);
    F = [F{:}]';

    function f = st2f(x)
        f = histc(x, edges);
    end
end
