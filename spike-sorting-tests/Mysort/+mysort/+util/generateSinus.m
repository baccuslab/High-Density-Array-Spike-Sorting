
function y = generateSinus(freq, phase, samples_per_second, samples)
L = samples;
Fs = samples_per_second;
T = 1/Fs;
t = (0:L-1)*T;                % Time vector

% 90Grad = pi/2;
phase = pi*phase/180;
% Sum of a 50 Hz sinusoid and a 120 Hz sinusoid
y = sin((2*pi*freq)*t + phase); 

