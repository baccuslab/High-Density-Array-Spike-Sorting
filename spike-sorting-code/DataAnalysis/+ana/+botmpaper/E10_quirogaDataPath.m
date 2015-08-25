function [dpath] = E10_quirogaDataPath()
    % Set here the path to the folder containing the simulated data of the 
    % Quiroga Benchmark usually found in the wave_clus package:
    % /wave-clus-2.0/Simulator/
    % The package can be downloaded here:
    % http://www.vis.caltech.edu/~rodri/Wave_clus/Wave_clus_home.htm
    
    pd = pdefs();
    if ~isempty(strfind(computer, 'WIN'))
%         dpath = 'C:\LocalData\WaveClus\Simulator\';
%         dpath = 'C:\LocalData\Quiroga\Simulator\';
        dpath = fullfile(pd.serverData, 'Quiroga', 'Simulator');
    else
        dpath = fullfile(pd.serverData, 'Quiroga', 'Simulator');
    end