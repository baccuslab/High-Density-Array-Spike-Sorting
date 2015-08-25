function f = interpolfun(y, method)
    % this function returns a functional for the waveform in
    % y that can interpolate at arbirary timepoints.
    %
    % sf = sincfun(y) returns a functional that interpolates the waveform
    % in y assuming that the values in yare equidistantly sampled with
    % unit distance between sample points
    %
    %
    % sf is a functional that takes one additional arguments: sf(t)
    % where t is the vector with requested sample points for the
    % interpolation
    % 
    % Example:
    %      sf = mysort.util.interpolfun(sin(0:.1:pi));
    %      figure;
    %      plot(sf(0:.01:pi));
    if nargin < 2
        method = 'pchip'; 
    end
    f = @(xi) interp1(1:length(y), y, xi, method, 0);    
end