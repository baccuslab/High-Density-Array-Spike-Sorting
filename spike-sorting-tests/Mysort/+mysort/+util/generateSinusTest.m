
Fs = 1024;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = 1024;                     % Length of signal
t = (0:L-1)*T;                % Time vector
% Sum of a 50 Hz sinusoid and a 120 Hz sinusoid
x = 1+ 0.7*sin(2*pi*50*t) + cos(2*pi*120*t) + cos(2*pi*505*t)+ cos(2*pi*512*t); 
y = x + 0*2*randn(size(t)); 

x2 = 1 + .7*mysort.util.generateSinus( 50,  0, Fs, L) ...
       +    mysort.util.generateSinus(120, 90, Fs, L) ...
       +    mysort.util.generateSinus(505, 90, Fs, L) ...
       +    mysort.util.generateSinus(512, 90, Fs, L);
   
norm(x-x2)

%%
x = mysort.util.generateSinus( 50,  0, Fs, L);
x2 = mysort.util.generateSinus( 50,  90, Fs, L);
norm(x-x2);

figure
plot(x, 'b')
hold on
plot(x2, 'r')
legend('sin 50Hz', 'cos 50Hz');
