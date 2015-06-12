
function [figs plotNames] = report(varargin)
    P.savepath = [];
    P.spikes = [];
    P.fetX = [];
    P.IDs = [];
    P.C = [];
    P.X = [];
    P.DH = [];
    P.SD = [];
    P.T = [];
    P.TIDs = [];
    P.nC = 1;
    P.gdf = [];
    P.srate = 32000;
    P.nameprefix = '';
    P.cutleft = 0;
    P = mysort.util.parseInputs(P, 'report', varargin);
    
    k = 1; plotNames = {};
    
    if isempty(P.DH) && ~isempty(P.X)
        P.DH = mysort.datafile.DataHandle(P.X, P.srate);
    end
    if ~isempty(P.DH)
        L = min(P.DH.Len, 100000);
    end
    
    if ~isempty(P.SD)
        plotNames{k} = 'spikeDetection';
        P.SD.plotDetection('start', 1, 'stopp', L);
        figs(k) = gcf; k=k+1;
    end
    
    if ~isempty(P.fetX) && ~isempty(P.IDs)
        plotNames{k} = 'clustering';
        mysort.plot.clustering(P.fetX, P.IDs); 
        figs(k) = gcf; k=k+1;
    end
    
    if ~isempty(P.spikes) && ~isempty(P.IDs) && ~isempty(P.T) && ~isempty(P.C)
        plotNames{k} = 'clusterProjections';
        mysort.plot.clusterProjection(P.spikes, P.IDs, P.T, P.C);
        figs(k) = gcf; k=k+1;
    end
    
    if ~isempty(P.spikes) && ~isempty(P.IDs) 
        plotNames{k} = 'spikes';
        mysort.plot.spikes(P.spikes, 'IDs', P.IDs, 'nC', P.nC);
        figs(k) = gcf; k=k+1;
    end
    
    if ~isempty(P.T) && ~isempty(P.IDs)    
        plotNames{k} = 'templates';
        mysort.plot.spikes(P.T, 'IDs', unique(P.IDs), 'nC', P.nC);
        figs(k) = gcf; k=k+1;

        plotNames{k} = 'templatesStacked';
        mysort.plot.spikes(P.T, 'IDs', unique(P.IDs), 'nC', P.nC, 'stacked', true);
        figs(k) = gcf; k=k+1;
    end
    
    if ~isempty(P.DH) && ~isempty(P.gdf) && ~isempty(P.T)
        plotNames{k} = 'sorting';
        mysort.plot.sorting(P.DH.getX('stopp', L), P.gdf(P.gdf(:,2)<L,:), P.T, 'template_IDs', P.TIDs, 'srate', P.DH.srate, 'cutleft', P.cutleft+1);
        figs(k) = gcf; k=k+1;
    end
    
    if ~isempty(P.spikes) && ~isempty(P.IDs) && ~isempty(P.C)
        plotNames{k} = 'intraClusterPCA';
        mysort.plot.clusters(P.spikes, P.IDs, P.C);
        figs(k) = gcf; k=k+1;
    end
     
    fprintf('Done.\n');
    
    % SAVE Plots
    if ~isempty(P.savepath)
        fprintf('Saving Plots... ');
        for k=1:length(figs)
            mysort.plot.savefig(figs(k), fullfile(P.savepath,  [P.nameprefix plotNames{k}]));
        end
    end