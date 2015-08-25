
function Stau = shiftMatrix(Tf, tau)
    n1 = Tf-abs(tau);
    Stau = diag(ones(1,n1), tau);
end