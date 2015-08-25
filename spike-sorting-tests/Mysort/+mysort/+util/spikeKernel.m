
function [k t spike2] = spikeKernel(spike1, spike2, iC, x0)
    % Computes the maximal xcorrelation between spike1 and spike2, the
    % corresponding shift and the shifted waveform
    % Input: 
    %   spike1, spike2 - Spike waveforms, channels as rows.
    %   iC             - inverse noise covariance function
    % Output:
    %   k       - maximal crosscorrelation
    %   t       - corresponding shift
    %   spike2  - shifted spike waveform of spike2 with shift t
    
    [nC L] = size(spike1);
    sRange = 1:L;   % DO NEVER CHANGE RESOLUTION AWAY FROM ONE SAMPLE
                    % THIS WILL CRIPPLE THE KERNEL, THE MAXIMUM OF THE
                    % AUTOCOVARIANCE FUNCTION IS NOT GARANTEED TO BE 
                    % 1 ANYMORE.
    if nargin < 3
        iC = [];
    end
    if nargin < 4
        x0 = 0; 
    end
    %assert(x1<x2, 'boundaries must be increasing!');
    import mysort.util.m2v
    import mysort.util.v2m

    spike2fun = mysort.util.sincfun(spike2);
    
    % DONT DO THIS! To find the maximum use the normal euclidean distance.
    % Use the Mahalanobis Distance after the shift is found, to calculate
    % the value of the crosscorrelation. Otherwise, the maximum of the
    % autocovariance function is not guaranteed to be at 0 !
%     if nargin < 3 || isempty(iC)
%         fun = @(t)   -sum(sum(spike1fun(sRange, 0).*spike2fun(sRange-t, 0)));
%     else
%         fun = @(t)   -m2v(spike1fun(sRange, 0))*iC*m2v(spike2fun(sRange-t, 0))';
%     end
    
    % Compute Shift firts, then energy:
    fun = @(t)   -sum(sum(spike1.*spike2fun(sRange-t, 0)));
    
    opt = [];
    opt.MaxFunEvals = 20;
    opt.TolX = 0.005;
    opt.Display = 'off';
    
    % Dont use fmibnd. It is much too likely that you get an obscure local
    % minimum. 
    %[t, k] = fminbnd(fun, x1, x2, opt);
    [t, k] = fminsearch(fun, x0, opt);
    
    if abs(t) < opt.TolX
        t = 0;
    end
    
    % Now calculate kernel
    if ~isempty(iC) || nargout ==3
        spike2 = spike2fun(sRange-t, 0);
    end
    if ~isempty(iC)
        k = m2v(spike1)*iC*m2v(spike2)';
    else
        k = -k; % switch sign, because fminsearch goes for minima!
    end
