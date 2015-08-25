
function [epochs keepIdx] = removeOverlapping(epochs, epochsRem)
    idx = zeros(size(epochs,1),1);
    count = 0; 
%     rmvIdx = [];
    for i=1:size(epochs,1)
        keep = 1;
        for k=1:size(epochsRem,1)
            if ((epochsRem(k,1) <= epochs(i,2)) && (epochsRem(k,2) >= epochs(i,1)))
                keep = 0;
                break
            end
        end
        if keep==1
            count = count +1;
            idx(count) = i;
%         elseif nargout>1
%             rmvIdx = [rmvIdx; i];
        end
    end
    keepIdx = idx(1:count);    
    epochs = epochs(keepIdx,:);
