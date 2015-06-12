function V = vChannelMul(V, n)
    warning('This function is depricated! Use mysort.wf.* instead!');
    % multiplies the channels of multichannel waveforms that are stored
    % in vector representation in V (single channel concatenated).
    % the number of channels is the number of elements in n.
    % each channel is multiplied by the repective number in n
    nC = length(n);
    nV = size(V,1);
    if size(n,2) > size(n,1); n=n'; end
    M = mysort.util.v2m(V, nC);
    M = M.*repmat(n, nV, size(M,2));
    V = mysort.util.m2v(M, nC);