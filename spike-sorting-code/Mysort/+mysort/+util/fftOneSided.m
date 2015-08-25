function [A freqs Y] = fftOneSided(x, Fs)
    L = length(x);                % Length of signal
    
    fNyq = Fs/2;
    
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    Y = fft(x,NFFT)/L;

    nFreqs = 1+length(Y)/2;
    A = 2*abs(Y(1:nFreqs));
    A([1 end]) = A([1 end])/2;
    freqs = linspace(0, fNyq, nFreqs);