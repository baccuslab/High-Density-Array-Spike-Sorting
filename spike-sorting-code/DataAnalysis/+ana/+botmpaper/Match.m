function R = Match(X,t, C, P)
if isempty(P)
    R = X*(t/C)';
else
    R = X*(P*(t/C)');
end
    
%     Eucl  = @(X, t)    sum(  (X-repmat(t, size(X,1), 1)).^2                                ,2);
%     Maha  = @(X, t, C) sum( ((X-repmat(t, size(X,1), 1))/C) .* (X-repmat(t, size(X,1), 1)) ,2);
%     Conv  = @(X, t)    X*t';
%     Match = @(X, t, C) X*(t/C)';
%     Botm  = @(X, t, C) X*(t/C)' - .5*t*(t/C)';