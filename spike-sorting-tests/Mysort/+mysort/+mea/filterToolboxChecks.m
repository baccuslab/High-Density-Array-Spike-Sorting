% This should run
hpf = 300;
lpf = 5000;
y.sr = 20000;
y.osf =1;
bp=fdesign.bandpass('n,f3dB1,f3dB2', 2, hpf, lpf, y.sr*y.osf);
%y.filters.bbp=butter(bp);
y.filters.bbp=design(bp,'butter','sosscalenorm','l1');


v = ver;
idx = find((strcmp('DSP System Toolbox', {v.Name})));
if isempty(idx)
    fprintf('DSP System Toolbox missing!');
else
    v(idx)
end

idx = find((strcmp('Signal Processing Toolbox', {v.Name})));
if isempty(idx)
    fprintf('Signal Processing Toolbox missing!');
else
    v(idx)
end

    