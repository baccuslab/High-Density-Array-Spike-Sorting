function sf = mInterpolfun(x,Y, method)
    % this function returns a functional for the multi channel waveform in
    % Y  (individual rows are waveforms) that can interpolate at arbirary
    % timepoints.
    %
    % sf = sincfun(Y) returns a functional that interpolates the waveform
    % in Y assuming that the values in Y are equidistantly sampled with
    % unit distance between sample points
    %
    % sf = sincfun(x,Y) uses the sample points in x associated with the y-
    % values in Y.
    %
    % sf is a functional that takes one additional arguments: sf(t)
    % where t is the vector with requested sample points for the
    % interpolation
    % 
    % Example:
    %      sf = mInterpolfun(sin(0:.1:pi));
    %      figure;
    %      plot(sf(0:.01:pi));
    if nargin < 3
        method = 'pchip'; 
        if nargin < 2
            Y=x;
            x=1:size(x,2);
        end
    end
    
    %assert(~any(diff(x) ~= 1), 'x values must be increasing by one!');
    
%     if size(y,2) == 1
%         y = y';
%     end
    
%     sf = @(xi) interp1(x,y,xi, method, 0);    
    sf = @(xi) interpHandle(xi);
    
    function Yi = interpHandle(xi)
        Yi = zeros(size(Y,1), length(xi));
        for i=1:size(Y,1)
            Yi(i,:) = interp1(x,Y(i,:), xi, method, 0);
        end
    end
end