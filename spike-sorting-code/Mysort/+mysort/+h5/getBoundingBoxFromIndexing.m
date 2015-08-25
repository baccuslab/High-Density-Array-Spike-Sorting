function [bb relIdx] = getBoundingBoxFromIndexing(maxLens, varargin)
    nDims = length(maxLens);
    nVarDims = length(varargin);
    assert(nVarDims <= nDims, 'cant index more dimensions then there are!')
    
    bb = zeros(2, nDims);
    relIdx = cell(1, nDims);
    bConsecutiveIndexing = 1;
    for i=1:nVarDims
        if ischar(varargin{i})
            assert(strcmp(varargin{i}, ':'), 'strange indexing!');
            bb(1,i) = 1;
            bb(2,i) = maxLens(i);
            relIdx{i} = varargin{i};
        else
            if any(diff(varargin{i}) ~= 1)
                bConsecutiveIndexing = 0;
                bb(1,i) = min(varargin{i});
                bb(2,i) = max(varargin{i});
                relIdx{i} = varargin{i} - bb(1,i) +1;
            else
                bb(1,i) = varargin{i}(1);
                bb(2,i) = varargin{i}(end);
                relIdx{i} = ':';
            end            
        end
    end
    for i=nVarDims+1:nDims
        bb(1,i) = 1;
        bb(2,i) = maxLens(i);
    end
    % if all indexing was consecutive, dont return the relative indexes
    % into the bounding box (that would be (:,:,...,:)
    if bConsecutiveIndexing
        relIdx = {};
    end