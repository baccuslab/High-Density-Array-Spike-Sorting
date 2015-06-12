
function merged = merge(epochs) 
    merged = [];
    if isempty(epochs)
        return
    end
    
    epochs = sortrows(epochs);
    
    merged(1,:) = epochs(1,:);
    k = 1;
    for i=2:size(epochs,1)
        if epochs(i,1) <= merged(k,2)
            merged(k,:) = [min(epochs(i,1), merged(k,1)) max(epochs(i,2), merged(k,2))];
        else
            k = k+1;
            merged(k,:) = epochs(i,:);
        end
    end