nC = 3;
s = 1:nC;
L = 5;
nT = 4;

vce = repmat(s, 1, L);
vce = repmat(vce, nT, 1) + repmat((0:nT-1)'*10, 1, nC*L);


vte = mysort.wf.vce2vte(vce, nC);

vce2 = mysort.wf.vte2vce(vte, nC);

assert(~any(vce(:)~=vce2(:)), 'error')

       