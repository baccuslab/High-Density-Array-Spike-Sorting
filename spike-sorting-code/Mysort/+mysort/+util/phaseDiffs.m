function phaseDiffs(x1, x2, Fs, varargin)
    assert(length(x1) == length(x2), 'x1 and x2 must be of same length');
    T = 1/Fs;                     % Sample time
    L = length(x1);                % Length of signal
    t = (0:L-1)*T;                % Time vector
    
    fNyq = Fs/2;
    if size(x1,1)<size(x1,2)
        X = [x1' x2'];
    else
        X = [x1 x2];
    end
    nFreqs = 1+length(X)/2;
    P = zeros(nFreqs, 2);
    for i=1:2
        x = X(:,i);
        NFFT = 2^nextpow2(L); % Next power of 2 from length of y
        Y = fft(x,NFFT)/L;
        f = 0:length(Y)-1;

        % Plot single-sided amplitude spectrum.
        
        YY = 2*abs(Y(1:nFreqs));
        YY([1 end]) = YY([1 end])/2;
        freqs = linspace(0, fNyq, nFreqs);
        phase = atan2(imag(Y(1:nFreqs)), real(Y(1:nFreqs)));
        phase(YY<10^-8) = nan;
        phase = 180*phase/pi;    
        P(:,i) = phase;
    end
    
%     ah(1) = subplot(2,1,1);
%     hold on
    d1 = abs(P(:,1)-P(:,2));
    d2 = abs(P(:,1)-P(:,2)+360);
    d3 = abs(P(:,1)-P(:,2)-360);
    plot(freqs, min([d1 d2 d3], [], 2), varargin{:}) 
    ylabel('phase diff(°)');
    xlabel('Frequency (Hz)')
    
%     linkaxes(ah, 'x');
% 