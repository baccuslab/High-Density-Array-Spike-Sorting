function config = configuration(varargin)
    config = mysort.SortingObject();
    config.name = 'default config';
    config.Tf = 70;
    config.srate = 20000;
    
    config.max_load_samples = 10*10^6;    
    config.max_simple_plot_channels = 16;    
    
    config = mysort.util.parseInputs(config, varargin, 'error');