%% WC
% Use wave_clus gui to sort the data. Only changes necessary to reproduce
% the 2004 results (table 2) are: 
% - use the set_parameters_simulation_default.m
% - reduce the cutting length of the waveforms to 
%   par.w_pre=9;                         % number of pre-event data points stored
%   par.w_post=25;                       % number of post-event data points stored
% - for Difficult 0.15 manual intervention is necessary to get 3 templates
% - for Difficult 0.2 two templates are correct, this was also in the 2004
%   paper
%
% Store all the clusters in the standard way
pd = pdefs();
simulator_path = ana.botmpaper.E10_quirogaDataPath();
if 0 
    addpath(genpath(fullfile(pd.serverData, 'Quiroga', 'wave_clus_2.0wb', 'Wave_clus')))
    wave_clus
    cd(simulator_path);
end

benchmarks = {'Easy1', 'Easy2', 'Difficult1', 'Difficult2'};
noiselevels = {'005', '01', '015', '02'};


%%
for b=1:length(benchmarks)
    for n=1:length(noiselevels)
        fprintf('\n###################\nProcessing %d - %d \n', b, n)
        quirogaFile = sprintf('C_%s_noise%s.mat', benchmarks{b}, noiselevels{n});
        data = load(fullfile(simulator_path, quirogaFile));
        gt_gdf = [data.spike_class{1}(:) data.spike_times{1}(:)];
        gt_gdf(:,2) = gt_gdf(:,2) + 30;
        
        wc_file = sprintf('times_C_%s_noise%s.mat', benchmarks{b}, noiselevels{n});
        sorted = load(fullfile(simulator_path, wc_file));
        sorted_gdf = sorted.cluster_class;
        sorted_gdf(sorted_gdf(:,1) == 0,:) = [];
        
        R = mysort.spiketrain.alignGDFs(gt_gdf, sorted_gdf, 10, 30, 10);
        mysort.plot.printEvaluationTable(R, 'ignoreEmptySorted', false);
        
                
    end
end
