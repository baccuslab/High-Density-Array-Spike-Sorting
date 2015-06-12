function f = filter_design(hpf, lpf, srate, forder, type)
    if nargin == 4
        type = 'butter';
    end
    if hpf == 0
        bp = fdesign.lowpass('n,f3dB', forder, lpf, srate);
    elseif lpf == srate/2
        bp = fdesign.highpass('n,f3dB', forder, hpf, srate);
    else
        bp = fdesign.bandpass('n,f3dB1,f3dB2', forder, hpf, lpf, srate);
    end
    
    f  = design(bp,type,'sosscalenorm','l1');