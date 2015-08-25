
function Tup = tResample(T, up, down, doNotCutTail)
    assert(int32(up) == up, 'up must be an integer!');
    assert(up>0, 'up must be bigger than 0!');
    assert(int32(down) == down, 'up must be an integer!');
    assert(down>0, 'up must be bigger than 0!');
    
    assert(ndims(T)<4 && ndims(T)>0, 'T must have 1, 2 or 3 dimensions!');
    
    if nargin < 4
        doNotCutTail = 0;
    end
    
    if up==1 && down == 1
        Tup = T;
        return
    end
    
    if ndims(T) < 3
        Tup = resample(T, up, down);
        if ~doNotCutTail
            Tup = Tup(1:end-(up-1));
        end
    else
        Tup = resample(T(:,:,1), up, down);
        Tup(end,size(T,2), size(T,3)) = 0;
        for i=2:size(T,3)
            Tup(:,:,i) = resample(T(:,:,i), up, down);
        end    
        if ~doNotCutTail
            % cut off interpolation behind the last real sample
            Tup = Tup(1:end-(up-1),:,:);
        end
    end    

