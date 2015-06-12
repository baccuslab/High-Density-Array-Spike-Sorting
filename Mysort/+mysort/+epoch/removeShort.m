
function epochs = removeShort(epochs, minLen)
    if isempty(epochs)
        return
    end

    len = epochs(:,2) - epochs(:,1) + 1;
    idx = zeros(size(len));
    for i=1:size(epochs,1)
        if len(i) >= minLen
            idx(i) = 1;
        end
    end
    epochs = epochs(idx>0,:);