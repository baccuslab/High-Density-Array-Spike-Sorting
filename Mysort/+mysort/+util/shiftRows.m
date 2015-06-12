function Y = shiftRows(X, tau, trunc)
    % shifts the rows of X by tau. truncates the sides of the new matrix
    % if trunc is supplied. truncates in a way, that the Y and X have the 
    % same dimensions. If a row has a shift of 0, and trunc is 1 it will
    % be unchanged.
    assert(nargin >= 2, 'X and tau need to be provided!');
    assert(nargin <= 3, 'Forth argument provided! Need shiftMCRows ?!');
    if nargin < 3; trunc = 0; end
    [rows, dim] = size(X);
    assert(rows == length(tau), 'For every row of X a shift value needs to be provided!');
    assert(trunc==0 || trunc==1, 'Trunc must be either true or false!');
    assert(size(tau,1) == 1 || size(tau,2) == 1, 'tau must not be a matrix!');
    minTau = min(tau);
    maxTau = max(tau);
    
    offset = max(0, -minTau);
    dimNew = dim + offset + max(0, maxTau);
    Y = zeros(rows, dimNew);
    
    for i=1:rows
        idxNew = offset+tau(i)+1:offset+tau(i)+dim;
        Y(i,idxNew) = X(i,:);
    end
    
    if trunc
        Y = Y(:, offset+1:offset+dim);
    end
    