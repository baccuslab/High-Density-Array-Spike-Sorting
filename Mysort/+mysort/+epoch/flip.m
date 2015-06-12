

function flipped = flip(epochs,len) 
% enough comment!!!

    flipped = [];
    if isempty(epochs)
        if nargin>1
            flipped = [1 len];
        end
        return        
    end

    if epochs(1,1) > 1
        flipped(1,:) = [1 epochs(1,1)-1];
    end
    flipped = [ flipped;
            [epochs(1:end-1,2)+1 epochs(2:end,1)-1] ];
    
    % Remove epochs of negative of zero length
    flipped(flipped(:,2)-flipped(:,1)<1,:) = [];
        
    if (nargin>1) && epochs(end,2)<len
        flipped = [flipped; 
            [epochs(end,2)+1 len] ];
    end    
    
