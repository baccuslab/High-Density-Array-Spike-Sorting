function b = filter_design_fir(hpf, lpf, srate, forder)
    if nargin < 4
        forder = 100;
    end
    nyq = srate/2;
    h = hpf/nyq;
    l = lpf/nyq;
    assert(h>=0 && h<=1, 'High pass is invalid!')
    assert(l>=0 && l<=1, 'Low pass is invalid!')
    assert(h < l, 'High pass must be smaller than low pass!')
    
    if h == 0
        bp = fdesign.lowpass('n,f3dB', forder, lpf, srate);
        error('lowpass is not yet implemented properly!')
    elseif lpf == srate/2
        bp = fdesign.highpass('n,f3dB', forder, hpf, srate);
        error('highpass is not yet implemented properly!')
    else
        amp   = [ 0     1      0];
        Upper = [ 0.001 1.01   0.001];
        Lower = [-0.001  .99  -0.001];
        Fc1 = h;
        Fc2 = l;
        b = fircls(forder, [0 Fc1 Fc2 1], amp, Upper, Lower);
    end
    
