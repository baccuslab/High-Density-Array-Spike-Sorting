classdef meanshiftSort < ana.botmpaper.harrisAnalysis.analysis
    properties
        PREC
    end
    
    methods
        %------------------------------------------------------------------
        function self = meanshiftSort(PREC)
            % Input argument PREC is optional but can be an
            % ana.botmpaper.harrisAnalysis.preprocess object
            self = self@ana.botmpaper.harrisAnalysis.analysis('MeanShift');
            if nargin == 1
                assert(isa(PREC, 'ana.botmpaper.harrisAnalysis.preprocess'), 'Wrong type of input argument!');
                self.PREC = PREC;
            else
                self.PREC = ana.botmpaper.harrisAnalysis.preprocess();
            end
        end
        %------------------------------------------------------------------
        function [fh figureNames] = makeFigures_(self, name, P)
            fh = []; p = 0; figureNames = {};
            R = self.getAllVariables(name);
            fh(p+1) = mysort.plot.figure('w', 1200, 'h', 1000); p=p+1;
                plot(R.S.spikePrewhitened.wfs')
                mysort.plot.figureTitle(['spikesPw - ' name]);
            figureNames{p} = 'prewhitenedSpikes';
            
            fh(p+1) = mysort.plot.figure('w', 1200, 'h', 1000); p=p+1;
                mysort.plot.clustering(R.S.spikeFeatures.X, R.S.clustering.ids, [], [], 'fh', fh(p), 'markerSize', 5)
                mysort.plot.figureTitle(['ms clustering - ' name]);
            figureNames{p} = 'msclustering';
                
            fh(p+1) = mysort.plot.figure('w', 1200, 'h', 1000); p=p+1;
                mysort.plot.clusterPCA(R.S.spikeFeatures.X, R.S.clustering.ids, [], 'fig', fh(p))
                mysort.plot.figureTitle(['ms clustering, PCA - ' name]);
            figureNames{p} = 'msclusteringPCA';
            
            % careful, these plots only work if all spikes were clustered
            fh(p+1) = mysort.plot.figure('w', 1200, 'h', 1000); p=p+1;
                mysort.plot.clustering(R.S.spikeFeatures.X, R.S.clusteringMerged.ids, [], [], 'fh', fh(p))
                mysort.plot.figureTitle(['clu merged - ' name]);
            figureNames{p} = 'mscluMerged';
                
            fh(p+1) = mysort.plot.figure('w', 1200, 'h', 1000); p=p+1;
                mysort.plot.clusterPCA(R.S.spikeFeatures.X, R.S.clusteringMerged.ids, [], 'fig', fh(p))
                mysort.plot.figureTitle(['clu merged, PCA - ' name]);
            figureNames{p} = 'mscluMergedPCA';         
            
            % careful, these plots only work if all spikes were clustered
            fh(p+1) = mysort.plot.figure('w', 1200, 'h', 1000); p=p+1;
                mysort.plot.clustering(R.S.spikeFeatures.X, R.S.clusteringMatched.ids, [], [], 'fh', fh(p))
                mysort.plot.figureTitle(['clu matched - ' name]);
            figureNames{p} = 'mscluMatched';
                
            fh(p+1) = mysort.plot.figure('w', 1200, 'h', 1000); p=p+1;
                mysort.plot.clusterPCA(R.S.spikeFeatures.X, R.S.clusteringMatched.ids, [], 'fig', fh(p))
                mysort.plot.figureTitle(['clu matched, PCA - ' name]);
            figureNames{p} = 'mscluMatchedPCA';   
        end
    end
    methods(Access=protected)
        %------------------------------------------------------------------
        function result = runTrial_(self, name) 
            Xfil = self.PREC.getVariable(name, 'Xfil');
            D = self.PREC.getVariable(name, 'D');
            
            DS = mysort.ds.Matrix(Xfil', D.srate, 'bla');
            dpath = fullfile(self.dataOutPrefix, 'MeanShiftSorting', D.name);

            P.spikeDetection.method = '-';
            P.spikeDetection.thr = 4.0;
            P.spikeDetection.minDist = 35;

            P.spikeCutting.Tf = 55;
            P.spikeCutting.cutLeft = 15;
            P.spikeCutting.maxSpikes = 200000;

            % Spike Alignment
            P.spikeAlignment.method = 'onUpsampledMean';
            P.spikeAlignment.maxSpikes = 200000;
            P.spikeAlignment.Tf = 25;
            P.spikeAlignment.cutLeft = 10;
            P.spikeAlignment.initAlignment = '-';
            P.spikeAlignment.maxIdx = P.spikeAlignment.cutLeft + 1;
            P.spikeAlignment.maxIterations = 30;

            % Noise estimation    
            P.noiseEstimation.minLength = 2000000;
            P.noiseEstimation.minDistFromSpikes = 60;

            % Feature extraction
            P.featureExtraction.Tf = 14;  
            P.featureExtraction.cutLeft = 7;
            P.featureExtraction.nDims = 6;

            % Clustering
            P.clustering.minSpikesPerCluster = 5;
            P.clustering.meanShiftBandWidthFactor = 1.2;
            P.clustering.maxSpikes = P.spikeAlignment.maxSpikes;

            P.mergeTemplates.atCorrelation = 0.99;
%             P.mergeTemplates.ifMaxDistSmaller = .5; % in noise std % this
            P.mergeTemplates.ifMaxRelDistSmallerPercent = .10;
            
            P.botm.cutLeft = P.spikeCutting.cutLeft;
            P.botm.run = 1;
            P.botm.Tf = P.spikeCutting.Tf;
            P.botm.prior = .000001;        

            [S P] = mysort.sorters.sort(DS, dpath, 'r2', P);
            detTs = S.spikeDetectionMerged.ts;
            gdf_clu = [S.clustering.ids detTs(S.spikeAligned.alignIdx)];
            ts = S.clusteringMatched.ts(:);
            gdf_match = [S.clusteringMatched.ids(:) ts];
            gdf_merged = [S.clusteringMerged.ids(:) ts];
            gdf_Botm = S.botm.gdf;
            
           
            result.S = S;
            result.P = P;      
            result.gdf_clu = gdf_clu;  
            result.gdf_match = gdf_match;
            result.gdf_merged = gdf_merged;
            result.gdf_Botm = gdf_Botm;
        end
        
    end
end