function R = Conv(X,t,P)
if isempty(P)
    R = X*t';
else
    R = X*(P*t');
end
    
%     Eucl  = @(X, t)    sum(  (X-repmat(t, size(X,1), 1)).^2                                ,2);
%     Maha  = @(X, t, C) sum( ((X-repmat(t, size(X,1), 1))/C) .* (X-repmat(t, size(X,1), 1)) ,2);
%     Conv  = @(X, t)    X*t';
%     Match = @(X, t, C) X*(t/C)';
%     Botm  = @(X, t, C) X*(t/C)' - .5*t*(t/C)';
