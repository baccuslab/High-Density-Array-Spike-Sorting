function Num = filter_gpu()
    p = mfilename('fullpath');
    [a b c] = fileparts(p);
    load(fullfile(a, 'gpuFilterCoeffcients2.mat'));
    
    %save('gpuFilterCoeffcients2.mat', 'Num');