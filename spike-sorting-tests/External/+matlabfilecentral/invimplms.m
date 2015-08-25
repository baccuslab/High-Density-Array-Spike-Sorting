function inv=invimplms(den,n,d)
% syntax inv=invimplms(den,n,d)
%         den - denominator impulse
%         n   - length of result
%         d   - delay of result
%         inv - inverse impulse response of length n with delay d
%
% Levinson-Durbin algorithm from Proakis and Manolokis p.865
%
% Author: Bob Cain, May 1, 2001 arcane[AT]arcanemethods[DOT]com

    m=xcorr(den,n-1);
    m=m(n:end);
    b=[den(d+1:-1:1);zeros(n-d-1,1)];
    inv=matlabfilecentral.toepsolve(m,b);





