function [T t1 t2] = cutTemplateToMiddleWithNDims(T, peakPosition, desiredDimsPerChannel)
    cutLeft = ceil(desiredDimsPerChannel/2)-1;
    cutRight = floor(desiredDimsPerChannel/2);
    
    t1 = peakPosition - cutLeft;
    t2 = peakPosition + cutRight;
    if t1 < 1
        t2 = t2 - t1 + 1;
        t1 = 1;
    elseif t2 > size(T,1)
        shift = t2 - size(T,1);
        t2 = size(T,1);
        t1 = t1 - shift;
    end
        
    T = T(t1:t2,:,:);
    assert(size(T,1) == desiredDimsPerChannel, 'Something went wrong!');