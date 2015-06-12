function sf = mSincfun(x,y)
    % this function returns a functional for the multi channel waveform in
    % x that can do an sinc interpolation at arbirary timepoints.
    %
    % sf = sincfun(y) returns a functional that interpolates the waveform
    % in y assuming that the values in y are equidistantly sampled with
    % unit distance between sample points
    %
    % sf = sincfun(x,y) uses the sample points in x associated with the y-
    % values in y.
    %
    % sf is a functional that takes one additional arguments: sf(t)
    % where t is the vector with requested sample points for the
    % interpolation
    % 
    % Example:
    %      sf = sincfun(sin(0:.1:pi));
    %      figure;
    %      plot(sf(0:.01:pi));
    if nargin == 1
        y=x;
        x=1:size(x,2);
    end
    
    %assert(~any(diff(x) ~= 1), 'x values must be increasing by one!');
    
    if size(y,2) == 1
        y = y';
    end
    
    sf = @(varargin) fun(x,y,varargin{:});
    
    function v = fun(x,y,t)
        % compute for every point in t(i) a sinc function that will be
        % multiplied with the values in y to yield the respective single
        % sample corresponding to t(i)
        S = mysort.util.sinc0(t(ones(size(x)),:) - x(ones(size(t)),:)');
        % since_ is a matrix containing the single interpolator sinc
        % functions as columns. They need to be of norm one
        nSinc_ = sqrt(sum(S.^2, 1));
        S = (S./nSinc_(ones(1, size(S,1)), :));
        v = y*S;
    end
end