
function M = t2m(T)
    warning('This function is depricated! Use mysort.wf.* instead!');
    [dim, channel, items] = size(T);
    M = zeros(items, dim*channel);
    for i=1:items
        tmp = squeeze(T(:,:,i));
        M(i,:) = tmp(:);
    end