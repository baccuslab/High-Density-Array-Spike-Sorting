function frequencyPhaseDiffsTest()
    Fs = 1000;                    % Sampling frequency
    T = 1/Fs;                     % Sample time
    L = 1024;                     % Length of signal
    t = (0:L-1)*T;                % Time vector
    fNyq= Fs/2;
    nFreqs = 1+L/2;
    freqs = linspace(0, fNyq, nFreqs);
        
    x1 = 0.7*sin(2*pi*30*freqs(2)*t) + sin(2*pi*60*freqs(2)*t); 
    x2 = 0.7*cos(2*pi*30*freqs(2)*(t+1)) + cos(2*pi*60*freqs(2)*t); 
    for i=1:200
        x1 = x1 + .5*cos(2*pi*(70+i*2)*freqs(2)*t + i/10);
    end
    for i=1:200
        x2 = x2 + .5*cos(2*pi*(70+i*2)*freqs(2)*t + i/5);
    end
    
%     y = x ;%+ 2*randn(size(t));     % Sinusoids plus noise
    

    figure;
    mysort.util.frequencyPhaseAnalysis(x1, Fs, '.-b');
    mysort.util.frequencyPhaseAnalysis(x2, Fs, '.-r');
    
    figure
    mysort.util.phaseDiffs(x1,x2,Fs, '.-');
    
