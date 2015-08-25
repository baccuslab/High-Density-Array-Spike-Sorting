
function v = m2v(m, nC)
    warning('This function is depricated! Use mysort.wf.* instead!');
    if nargin == 1
        % opposite of v2m: x := v2m(m2v(x),size(x,1))
        v = m';
        v = v(:)';
    else
        nVecs = size(m,1)/nC;
        nConcatDims = size(m,2)*nC;
        v = zeros(nVecs, nConcatDims);
        for i=1:nVecs
            idx = (i-1)*nC+1:i*nC;
            v(i,:) = mysort.wf.m2v(m(idx,:));
        end
    end
end