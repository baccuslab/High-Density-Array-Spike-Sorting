
function D = maha(X,Y, SIG)
    D = zeros(size(X,1),1);
    IS = inv(SIG);
    TD = X-Y;
    for i=1:size(X,1)
        D(i) = sqrt(TD(i,:) * IS * TD(i,:)');
    end