function hinv=toepsolve(r,q)
% Solve Toeplitz system of equations.
%    Solves R*hinv = q, where R is the symmetric Toeplitz matrix
%    whos first column is r
%    Assumes all inputs are real
%    Inputs:  
%       r - first column of Toeplitz matrix, length n
%       q - rhs vector, length n
%    Outputs:
%       hinv - length n solution
%
%   Algorithm from Roberts & Mullis, p.233
%
%   Author: T. Krauss, Sept 10, 1997
%
%   Modified: R. Cain, Dec 16, 2004 to remove a pair of transposes
%             that caused errors.

    n=length(q);
    a=zeros(n+1,2);
    a(1,1)=1;

    hinv=zeros(n,1);
    hinv(1)=q(1)/r(1);

    alpha=r(1);
    c=1;
    d=2;

    for k=1:n-1,
       a(k+1,c)=0;
       a(1,d)=1;
       beta=0;
       j=1:k;
       beta=sum(r(k+2-j).*a(j,c))/alpha;
       a(j+1,d)=a(j+1,c)-beta*a(k+1-j,c);
       alpha=alpha*(1-beta^2);
       hinv(k+1,1)=(q(k+1)-sum(r(k+2-j).*hinv(j,1)))/alpha;
       hinv(j)=hinv(j)+a(k+2-j,d)*hinv(k+1);
       temp=c;
       c=d;
       d=temp;
    end