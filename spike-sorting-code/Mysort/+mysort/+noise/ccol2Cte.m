function Cte = ccol2Cte(Ccol, Tf_aim)
    nC = size(Ccol,2);
    Tf = size(Ccol,1)/nC;
%     maxLag = Tf-1;
    if nargin == 1
        Tf_aim = Tf;
    end
    
    if Tf_aim < Tf
        Ccol = Ccol(1:(Tf_aim*nC),:);
    elseif Tf_aim == Tf
        % do nothing
    else
        error('not implemented! Dont do this, zero padding gives bad results!');
    end
    
    xc = mysort.noise.ccol2xcorr(Ccol);
    Cte = mysort.noise.xcorr2Cte(xc);
    
    