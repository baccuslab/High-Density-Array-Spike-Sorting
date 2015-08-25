s = 1:5;
nC = 3;
nT = 4;

vte = repmat(s, 1, nC);
vte = repmat(vte, nT, 1) + repmat((0:nT-1)'*10, 1, nC*length(s));


vce = mysort.wf.vte2vce(vte, nC);
vce2 = mysort.util.embedTime2embedChan(vte, nC);
vte_ = mysort.wf.vce2vte(vce, nC);
vte2_ = mysort.wf.vce2vte(vce2, nC);

vte
vte_
vte2_
vce
vce2