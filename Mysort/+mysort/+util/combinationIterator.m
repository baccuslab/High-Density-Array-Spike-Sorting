function varargout = combinationIterator(varargin)
    cardinalities = cellfun(@(x) length(x), varargin);
    allCombinations = mysort.util.buildCombinatorics(cardinalities);
    nP = length(cardinalities);
    varargout = cell(1, nP);
    for k=1:nP
        varargout{k} = [varargin{k}(allCombinations(:,k))];
    end