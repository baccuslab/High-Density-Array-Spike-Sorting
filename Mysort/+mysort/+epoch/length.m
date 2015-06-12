function L = length(epochs)
    if isempty(epochs)
        L = [];
        return
    end
    L = epochs(:,2)-epochs(:,1)+1;