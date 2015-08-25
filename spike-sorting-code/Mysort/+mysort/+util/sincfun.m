function sf = sincfun(x,y)
    error('this function is depricated. use mysort.wf.* instead!');
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
    % sf is a functional that takes two additional arguments: sf(t,d)
    % where t is the vector with requested sample points for the
    % interpolation and d is the derivative. d=0 means no derivation.
    % 
    % Example:
    %      sf = sincfun(sin(0:.1:pi));
    %      figure;
    %      plot(sf(0:.01:pi, 0));
    if nargin == 1
        y=x;
        x=1:size(x,2);
    end
    
    %assert(~any(diff(x) ~= 1), 'x values must be increasing by one!');
    
    if size(y,2) == 1
        y = y';
    end
    
    sf = @(varargin) fun(x,y,varargin{:});
    
    function v = fun(x,y,t,varargin)
        nC = size(y,1);
        L  = length(t);
        v = zeros(nC, L);
        sinc_ = mysort.util.sinc(t(ones(size(x)),:) - x(ones(size(t)),:)',varargin{:});
        sinc_ = sinc_/norm(sinc_);
        for c=1:nC
            v(c, :) = (-1)^d * y(c,:)*sinc_;
        end
    end
end