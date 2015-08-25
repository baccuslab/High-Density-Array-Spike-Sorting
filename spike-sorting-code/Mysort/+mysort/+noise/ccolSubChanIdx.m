function ccol = ccolSubChanIdx(CCol, chanidx, maxLag_req)
nC = size(CCol,2);
Tf = size(CCol,1)/nC;
maxLag = Tf-1;

if nargin < 3
    maxLag_req = maxLag;
end

nC_req = length(chanidx);
xidx = zeros(Tf*nC_req, 1);

for i=1:Tf
    xidx((i-1)*nC_req+1:i*nC_req, 1) = chanidx + (i-1)*nC;
end
ccol = CCol(xidx, chanidx);

if maxLag_req < maxLag
    ccol = ccol(1:nC_req*(maxLag_req+1),:);
end