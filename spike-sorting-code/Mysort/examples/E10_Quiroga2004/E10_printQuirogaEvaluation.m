function E10_printQuirogaEvaluation(benchmark)
    Q = E10_QuirogaPerformanceIn2004Paper();
    
    b = benchmark2Idx(benchmark);
    
    fprintf('\nEvaluation reported in 2004 publication Quiroga et al:\n');
    mysort.plot.printTable(Q{b}.Table, 'colLabel', Q{b}.colLabel,...
                                       'rowLabel', Q{b}.rowLabel,...
                           'topLeftLabel', Q{b}.name, 'printColSum', 1);
        
    %----------------------------------------------------------------------
    function b = benchmark2Idx(benchmark)
        if ischar(benchmark)
            if strcmp(benchmark, 'Easy 1') || strcmp(benchmark, 'Easy1')
                b = 1;
            elseif strcmp(benchmark, 'Easy 2') || strcmp(benchmark, 'Easy2')
                b = 2;
            elseif strcmp(benchmark, 'Difficult 1') || strcmp(benchmark, 'Difficult1')
                b = 3;
            else
                b = 4;
            end
        else
            b = benchmark;
        end
    end
end