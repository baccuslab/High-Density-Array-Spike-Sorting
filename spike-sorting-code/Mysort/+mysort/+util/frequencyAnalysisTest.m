function frequencyAnalysisTest()
    Fs = 1000;                    % Sampling frequency
    T = 1/Fs;                     % Sample time
    L = 2*1024;                     % Length of signal
    t = (0:L-1)*T;                % Time vector
    fNyq= Fs/2;
    nFreqs = 1+L/2;
    freqs = linspace(0, fNyq, nFreqs);
        
    x = 0.7*sin(2*pi*50*freqs(2)*t) + sin(2*pi*120*freqs(2)*t); 
    y = x ;%+ 2*randn(size(t));     % Sinusoids plus noise
    figure;
    mysort.util.frequencyAnalysis(y, Fs, 'onesided', '.-');
    figure;
    mysort.util.frequencyAnalysis(y, Fs, 'twosided', '.-');