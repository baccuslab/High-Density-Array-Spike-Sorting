
function Tup = resampleTensor(T, up, down)
    warning('Dont use this function (mysort.util.resampleTensor), use mysort.wf.tResample instead!');
    assert(int32(up) == up, 'up must be an integer!');
    assert(up>0, 'up must be bigger than 0!');
    assert(int32(down) == down, 'up must be an integer!');
    assert(down>0, 'up must be bigger than 0!');
    
    assert(ndims(T)<4 && ndims(T)>0, 'T must have 1, 2 or 3 dimensions!');
    if ndims(T) < 3
        Tup = resample(T, up, down);
        Tup = Tup(1:end-(up-1));
    else
        Tup = resample(T(:,:,1), up, down);
        Tup(end,size(T,2), size(T,3)) = 0;
        for i=2:size(T,3)
            Tup(:,:,i) = resample(T(:,:,i), up, down);
        end    
        % cut off interpolation behind the last real sample
        Tup = Tup(1:end-(up-1),:,:);
    end    

