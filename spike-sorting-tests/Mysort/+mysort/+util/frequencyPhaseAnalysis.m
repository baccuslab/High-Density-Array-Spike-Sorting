
function [A freqs uw_phase phase_delay] = frequencyPhaseAnalysis(x, Fs, amplitude_cutoff, varargin)
% claculates the one or two sided frequency spectrum of data in x
% inputs:
%       x - data, single channel
%      Fs - sampling rate in samples per second
%     all other arguments will be handed to the plot function
% example:
%   x = randn(10000,1);
%   y = filter([1 1 1]/3, 1, x);
    if nargin < 3 || isempty(amplitude_cutoff)
        amplitude_cutoff = 10^-8;
    end

    [A freqs Y] = mysort.util.fftOneSided(x, Fs);

    [phase uw_phase phase_delay] = mysort.util.phaseAnalysis(Y, A, freqs, amplitude_cutoff);
    
    N = 3;
    ah(1) = subplot(N,1,1);
    hold on
    plot(freqs, A, varargin{:}) 
    title('One-Sided Amplitude Spectrum of x(t)')
    ylabel('|Y(f)|')
    
    
    ah(2) = subplot(N,1,2);
    hold on

    plot(freqs, uw_phase, varargin{:});
%     set(ah(2), 'ylim', [-180 180]);
    ylabel('phase (°)');
    xlabel('Frequency (Hz)')
    title('Phase for non-zero frequencies');
    
    ah(3) = subplot(N,1,3);
    hold on
    plot(freqs, phase_delay, varargin{:});    
    ylabel('phase delay');
    linkaxes(ah, 'x');
