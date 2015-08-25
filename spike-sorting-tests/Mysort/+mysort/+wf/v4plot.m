function [V rangeidx] = v4plot(V, nC)
    V = mysort.wf.m2v([mysort.wf.v2m(V, nC) nan(size(V,1)*nC,1)], nC);
    if nargout == 2
        rangeidx = mysort.wf.vSubChannelIdx(V, nC);
    end