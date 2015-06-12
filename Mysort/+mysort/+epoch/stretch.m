
function epochs = stretch(epochs, len)
    if isempty(epochs)
        return
    end
    epochs(:,1) = epochs(:,1)-len;
    epochs(:,2) = epochs(:,2)+len;