
function [YY freqs PHASE] = frequencyAnalysis(x, Fs, mode, varargin)
% claculates the one or two sided frequency spectrum of data in x
% inputs:
%       x - data, single channel
%      Fs - sampling rate in samples per second
%    mode - either 'onesided' or 'twosided'
%     all other arguments will be handed to the plot function
% example:
%   x = randn(10000,1);
%   y = filter([1 1 1]/3, 1, x);
%   figure
%   mysort.util.frequencyAnalysis(x, 10000, 'onesided', '.-b')
%   hold on
%   mysort.util.frequencyAnalysis(y, 10000, 'onesided', 'r')
    if nargin < 3 || isempty(mode);
        mode = 'onesided';
    end
    T = 1/Fs;                     % Sample time
    L = length(x);                % Length of signal
    t = (0:L-1)*T;                % Time vector
    
    fNyq = Fs/2;
    
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    Y = fft(x,NFFT)/L;
    f = 0:length(Y)-1;

    if strcmp(mode,'twosided')
        % Plot single-sided amplitude spectrum.
        plot(f,2*abs(Y), varargin{:}) 
        PHASE = imag(Y);
        title('Two-Sided Amplitude Spectrum of y(t)')
        xlabel('Frequency (Hz)')
        ylabel('|Y(f)|')
    elseif strcmp(mode,'onesided')
        % Plot single-sided amplitude spectrum.
        nFreqs = 1+length(Y)/2;
        YY = 2*abs(Y(1:nFreqs));
        PHASE = imag(Y(1:nFreqs));
        YY([1 end]) = YY([1 end])/2;
        freqs = linspace(0, fNyq, nFreqs);
        plot(freqs, YY, varargin{:}) 
        title('One-Sided Amplitude Spectrum of x(t)')
        xlabel('Frequency (Hz)')
        ylabel('|Y(f)|')
    else
        error('mode unknown: %s', mode);
    end