function [out occurences idx_occurences] = findCommonClusters(II, varargin)
% Input vector II must be a multi column vector with clustering indices
% for in each column.
% The length of each colum corresponds to the number of clustered units.
assert(size(II, 2) > 1, 'Input vector must have at least two columns!');

P.doPlot = false;
P = mysort.util.parseInputs(P, varargin, 'error');

[Uq Nq Gq] = unique(II, 'rows'); % Uq: unique rows, Uq = II(Nq,:), II = Uq(Gq,:)
[occurences_ x] = hist(Gq, length(Uq));
[occurences sort_idx] = sort(occurences_,'descend');



out = {};
for i = 1:length(sort_idx)
    out{i} = find(Gq == sort_idx(i));
end
idx_occurences = Uq(sort_idx,:);

if P.doPlot
    figure; scatter(idx_occurences(:,1), idx_occurences(:,2), 200*occurences)
    text(idx_occurences(:,1),idx_occurences(:,2), num2str(occurences') );
end

end