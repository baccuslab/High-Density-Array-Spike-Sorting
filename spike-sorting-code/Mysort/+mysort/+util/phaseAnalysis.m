function [phase uw_phase phase_delay] = phaseAnalysis(Y, YY, freqs, amplitude_cutoff)
    if nargin < 3
        amplitude_cutoff = 10^-8;
    end
    nFreqs = length(freqs);
    phase = atan2(imag(Y(1:nFreqs)), real(Y(1:nFreqs)));
    phase(YY<amplitude_cutoff) = nan;
    uw_phase = unwrap(phase, 2);
    uw_phase = 180*uw_phase/pi;

    phase_delay = nan(size(uw_phase));
    for i=2:length(uw_phase)
        phase_delay(i) = -uw_phase(i)/freqs(i);
    end