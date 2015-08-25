% -------------------------------------------------------------------------
% Corrects rates of a multivariate Poisson distribution for a constant
% refractory period of 5 ms.
% lambda - Mean rates of independent Poisson variables
% A - Mix matrix
% -------------------------------------------------------------------------
function lambda = mvpoissrefractoryrates(lambda, A)

lambda_orig = lambda;
refcost = @(x)mse(x, A, lambda_orig);
lambda = fminsearch(refcost, lambda);

if (any(lambda<0))
    [v, vindex] = min((lambda<0) .* lambda);
    error(['mvpoissrefractoryrates: Cannot generate requested covariances. Negative variance: ', num2str(v), ' at index ', int2str(vindex), '.']);
end

% Calculates the mean square error caused by the refractory period
function val = mse(x, A, lambda)

rates = A * x;

% Square polynom that interpolates the reduction in rate caused by the
% refractory period
refpoly = [-0.0044 0.9945 0.0084];
rates_ref = polyval(refpoly, rates);

% Calculate resulting effect on covariances
x_ref = zeros(size(x));
for i = 1:length(x_ref)
    indices = A(:, i) > 0;
    x_ref(i) = x(i) .* min(rates_ref(indices) ./ rates(indices));
end

val = mean((x_ref - lambda) .^ 2) - sum(x_ref(x_ref < 0));
