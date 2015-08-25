function bbp = filter_design_fpga(hpf, lpf, srate, forder)
    % Designs an IIR filter df2, that is used on the FPGA
    
    if nargin < 4
        forder = 2;
    end
    if nargin < 3
        srate = 20000;
    end
    if nargin < 2
        lpf = 3000;
    end
    if nargin < 1
        hpf = 500;
    end
    
    bp  = fdesign.bandpass('n,f3dB1,f3dB2', forder, hpf, lpf, srate);
    bbp = butter(bp);
    
