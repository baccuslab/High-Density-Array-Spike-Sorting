
function A = buildCombinatorics(lens, changeLastColumnFastest)
%     nCols = length(lens);
%     nRows = prod(lens);
%     C = zeros(nRows, nCols);
%  
%     for i=0:nRows-1
%         C(i+1,:) = mysort.util.getNumberInStrangeBasis(i, lens);
%     end
%     C = C+1;


    % Alternative that is faster !!

    if iscell(lens)
        % use the actual items in each cell in lens to build combinatorics
        if nargin == 1 || ~changeLastColumnFastest
            [X{1:length(lens)}] = ndgrid(lens{:});
            A = reshape(cat(length(lens)+1,X{:}),[],length(lens)); 
        else
            lens = lens(end:-1:1);
            [X{1:length(lens)}] = ndgrid(lens{:});
            A = reshape(cat(length(lens)+1,X{:}),[],length(lens));
            A = fliplr(A);
        end
    else
        lens = lens(:)';
        
        if nargin == 1 || ~changeLastColumnFastest
            A = helper(lens);
        else
            lens = fliplr(lens);
            A = helper(lens);
            A = fliplr(A);
        end
    end
    
    function X = helper(L)
        lensAsArrays = arrayfun(@(x) 1:x, L, 'uniformOutput', false);
        [X{1:length(L)}] = ndgrid(lensAsArrays{:});
        X = reshape(cat(length(L)+1,X{:}),[],length(L));        
    end
end