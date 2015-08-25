classdef clusterEval < ana.botmpaper.harrisAnalysis.analysis
    properties
        PREC
        PRECGT
        MS
        
    end
    
    methods
        %------------------------------------------------------------------
        function self = clusterEval()
            self = self@ana.botmpaper.harrisAnalysis.analysis('ClusterEval');
            self.PREC = ana.botmpaper.harrisAnalysis.preprocess();
            self.PRECGT = ana.botmpaper.harrisAnalysis.preprocessGroundTruth();
            self.MS = ana.botmpaper.harrisAnalysis.meanshiftSort();
        end
        %------------------------------------------------------------------
        function [fh figureNames] = makeFigures_(self, name, P)
            fh = []; p = 0; figureNames = {};
            R = self.getAllVariables(name);
            Rms  = self.MS.getAllVariables(name);
            names = {};
            SpikeTypes = zeros(length(R.E), 4);
            for i=1:length(R.E)
                names{i} = R.E{i}.sortingName;
                [fh_ figureNames_] = self.makeEvalPlots(R.E{i}, name, Rms);
                SpikeTypes(i,:) = [size(R.E{i}.spikesETP,1) size(R.E{i}.spikesEFP,1) size(R.E{i}.spikesEFN,1) size(R.E{i}.spikesECL,1)];
                fh = [fh fh_];
                figureNames = [figureNames figureNames_];
            end
            p = length(fh);
            
            fh(p+1) = mysort.plot.figure('w', 1200, 'h', 1000); p=p+1;
              nR=2; nC = 3; ap=0; ah =[];
              ap =ap+1; ah(ap) = subplot(nR,nC,ap); 
              h = bar(SpikeTypes(:,1));
              set(h, 'facecolor', 'g')
              set(ah(ap), 'xticklabel', names)
              title('True Positives');
              
              ap =ap+1; ah(ap) = subplot(nR,nC,ap); 
              h = bar(SpikeTypes(:,2));
              set(h, 'facecolor', 'b')
              set(ah(ap), 'xticklabel', names)
              title('False Positives');
              
              ap =ap+1; ah(ap) = subplot(nR,nC,ap); 

              ap =ap+1; ah(ap) = subplot(nR,nC,ap); 
              h = bar(SpikeTypes(:,4));
              set(h, 'facecolor', 'r')
              set(ah(ap), 'xticklabel', names)
              title('Clustering Errors');
              
              ap =ap+1; ah(ap) = subplot(nR,nC,ap); 
              h = bar(SpikeTypes(:,3));
              set(h, 'facecolor', 'k')
              set(ah(ap), 'xticklabel', names)
              title('False Negatives');           
              
              ap =ap+1; ah(ap) = subplot(nR,nC,ap); 
              h = bar(sum(SpikeTypes(:,2:end),2));
              set(h, 'facecolor', [.5 .5 .5])
              set(ah(ap), 'xticklabel', names)
              title('All Errors'); 
              mysort.plot.figureTitle([name ' - Error Overview']);
            figureNames{p} = 'ErrorTypesOverview';
%             
% 
%             fh(p+1) = mysort.plot.figure('w', 1200, 'h', 1000); p=p+1;
%                
%             figureNames{p} = '';   
                
        end
 %------------------------------------------------------------------
        function [fh figNames] = makeEvalPlots(self, E, name, Rms)
            p=0; ap = 0;
            fh(p+1) = mysort.plot.figure('w', 1200, 'h', 1000); p=p+1;
                ah(ap+1) = subplot(5,2,1); ap =ap+1;
                plot(E.spikesIFP');
                title(sprintf('false positives (%d)',size(E.spikesIFP,1)))          
                ah(ap+1) = subplot(5,2,3); ap =ap+1;
                plot(E.spikesICL');
                title(sprintf('clustering errors (%d)',size(E.spikesICL,1)))          
                ah(ap+1) = subplot(5,2,5); ap =ap+1;
                plot(E.spikesIFN');
                title(sprintf('false negatives (%d)',size(E.spikesIFN,1)))          
                ah(ap+1) = subplot(5,2,7); ap =ap+1;
                plot(E.spikesITP');
                title(sprintf('true positives (%d)',size(E.spikesITP,1)))                
                ah(ap+1) = subplot(5,2,9); ap =ap+1;
                plot(median(E.spikesITP,1), 'g', 'linewidth', 2); hold on
                plot(median(E.spikesIFP,1), 'b', 'linewidth', 2);
                plot(median(E.spikesICL,1), 'r', 'linewidth', 2);
                plot(median(E.spikesIFN,1), 'k', 'linewidth', 2);                
                
                ah(ap+1) = subplot(5,2,2); ap =ap+1;
                plot(E.spikesEFP');
                title(sprintf('false positives (%d)',size(E.spikesEFP,1)))          
                ah(ap+1) = subplot(5,2,4); ap =ap+1;
                plot(E.spikesECL');
                title(sprintf('clustering errors (%d)',size(E.spikesECL,1)))          
                ah(ap+1) = subplot(5,2,6); ap =ap+1;
                plot(E.spikesEFN');
                title(sprintf('false negatives (%d)',size(E.spikesEFN,1)))  
                ah(ap+1) = subplot(5,2,8); ap =ap+1;
                plot(E.spikesETP');
                title(sprintf('true positives (%d)',size(E.spikesETP,1)))                                
                ah(ap+1) = subplot(5,2,10); ap =ap+1;
                plot(median(E.spikesETP,1), 'g', 'linewidth', 3);hold on
                plot(median(E.spikesEFP,1), 'b', 'linewidth', 2); 
                plot(median(E.spikesECL,1), 'r', 'linewidth', 2);
                plot(median(E.spikesEFN,1), 'k', 'linewidth', 2);
                                
                axis(ah, 'tight')
                legend('TP', 'FP', 'CL', 'FN', 'location', 'northeast')
                title('STAs')       
                mysort.plot.figureTitle([name ' - ' E.sortingName ' - Error Waveforms']);
            figNames{p} = [E.sortingName 'ErrorWaveforms'];
            
            fh(p+1) = mysort.plot.figure('w', 1400, 'h', 1000); p=p+1;
                cluP = mysort.plot.clustering(E.spikeFeatures, E.sgdf(:,1), [], [], 'errorTypes', E.eT, 'fh', fh(p));  
                mysort.plot.legend(cluP.axesHandles(2), {{'o', 'color', [1 0 0], 'markersize', 7, 'linewidth', 2},...
                                         {'o', 'color', [0 0 1], 'markersize', 7, 'linewidth', 2},...
                                         {'.', 'color', mysort.plot.vectorColor(E.targetUnitID)}},...
                                        {'CL/FN', 'FP', 'targetNeuron'});
                mysort.plot.figureTitle([name ' - ' E.sortingName ' - Errors In Clusterspace']);
                
            figNames{p} = [E.sortingName 'ErrorsInClusterspace'];
           
            fh(p+1) = mysort.plot.figure('w', 1400, 'h', 1000); p=p+1;
                mysort.plot.clusterPCA(E.spikeFeatures, E.sgdf(:,1), [], 'errorTypes', E.eT, 'fig', fh(p))  
                mysort.plot.figureTitle([name ' - ' E.sortingName ' - Errors In ClusterPCA']);
            figNames{p} = [E.sortingName 'ErrorsInClusterPCA'];                        
%             mysort.plot.clusterProjection(S.spikePrewhitened.wfs(u3idx,:), u3Et, [], [], 'plotTitleDiff', 1, 'plotTitleRsq', 1) 
%             mysort.plot.clusterProjection(S.spikeFeatures.X(u3idx,:), u3Et, [], [], 'plotTitleDiff', 1, 'plotTitleRsq', 1)             
        end        
        %------------------------------------------------------------------
        function E = makeEval(self, tgdf, sgdf, I, X)
            def = mysort.util.defs();
            
%             round([tgdf(1:100,2) sgdf(1:100,:)])
            E.sgdf = sgdf;
            E.tgdf = tgdf;
            E.Rinv = mysort.spiketrain.alignGDFs(sgdf, tgdf, 10, 10, 10);
            E.rstrinv = mysort.plot.printEvaluationTable(E.Rinv, 'ignoreEmptySorted', false); disp(E.rstrinv)     

            E.R = mysort.spiketrain.alignGDFs(tgdf, sgdf, 10, 10, 10);
            E.rstr = mysort.plot.printEvaluationTable(E.R, 'ignoreEmptySorted', false); disp(E.rstr);
            
    %         mysort.util.logToFile(fullfile(dpath, 'rstr.txt'), ['\n' rstr]);        
    
            targetIdx = E.R.k2f(1);
            uIDs = unique(sgdf(:,1));
            E.targetUnitID = targetIdx; 
            if uIDs(1) == 0
                % If the IDs start with 0, the index is shifted by one
                E.targetUnitID = E.targetUnitID-1;
            end
            
            
            E.eT = mysort.spiketrain.alignGetErrorTypes(E.R, sgdf, E.targetUnitID, uIDs);

            E.targetUnitSpikeIdx = find(sgdf(:,1)==E.targetUnitID);
            E.targetUnitEt = E.R.spikeLabel2{targetIdx};
            
            E.targetUnitTPIdx = E.targetUnitSpikeIdx(E.targetUnitEt==def.tp);
            E.TPts = sgdf(E.targetUnitTPIdx, 2);
            
            E.targetUnitFPIdx = E.targetUnitSpikeIdx(E.targetUnitEt==def.fp);
            E.FPts = sgdf(E.targetUnitFPIdx, 2);

            E.targetUnitCLIdx = find(E.R.spikeLabel1{1}==def.cl);
            E.CLts = E.tgdf(E.targetUnitCLIdx, 2);
            
            E.targetUnitFNIdx = find(E.R.spikeLabel1{1}==def.fn);
            E.FNts = E.tgdf(E.targetUnitFNIdx, 2);
            
            E.spikesITP = mysort.epoch.extractWaveform(I, round([E.TPts-10 E.TPts+40]));
            E.spikesIFP = mysort.epoch.extractWaveform(I, round([E.FPts-10 E.FPts+40]));
            E.spikesICL = mysort.epoch.extractWaveform(I, round([E.CLts-10 E.CLts+40]));
            E.spikesIFN = mysort.epoch.extractWaveform(I, round([E.FNts-10 E.FNts+40]));
            
            E.spikesETP = mysort.epoch.extractWaveform(X, round([E.TPts-10 E.TPts+40]));
            E.spikesEFP = mysort.epoch.extractWaveform(X, round([E.FPts-10 E.FPts+40]));
            E.spikesECL = mysort.epoch.extractWaveform(X, round([E.CLts-10 E.CLts+40]));
            E.spikesEFN = mysort.epoch.extractWaveform(X, round([E.FNts-10 E.FNts+40]));
        end 
    end
    methods(Access=protected)
        %------------------------------------------------------------------
        function R = runTrial_(self, name)
            Rgt = self.PRECGT.getAllVariables(name);
            if isempty(Rgt) 
                R = 'No ground truth available';
                return
            end
            Rms  = self.MS.getAllVariables(name);
            Rpre = self.PREC.getAllVariables(name);
            
            
            detTs = Rms.S.spikeDetectionMerged.ts;
            gdf_ms = [Rms.S.clustering.ids detTs(Rms.S.spikeAligned.alignIdx)];
            gdf_match = [Rms.S.clusteringMatched.ids(:) Rms.S.clusteringMatched.ts(:)];
            gdf_botm = Rms.S.botm.gdf;
            gdf_botm(:,2) = gdf_botm(:,2) + floor(Rms.P.botm.Tf/2)-10;
            sortings = {gdf_ms gdf_match gdf_botm};
            sortingNames = {'MSClustering' 'MSMatched' 'BOTM'};
            E = cell(length(sortings),1);
            for i=1:length(sortings)
                E{i,1} = self.makeEval(Rgt.tgdfE, sortings{i}, Rpre.I, Rpre.Xfil);
                E{i,1}.sortingName = sortingNames{i};
                if strcmp(sortingNames{i}, 'BOTM')
                    cl = Rms.P.featureExtraction.cutLeft;
                    Tf = Rms.P.featureExtraction.Tf;
                    spikes = mysort.epoch.extractWaveform(Rpre.Xfil, [gdf_botm(:,2)-cl gdf_botm(:,2)+Tf-cl-1]);
                    spikes = spikes/Rms.S.spikePrewhitened.U;
                    E{i,1}.spikeFeatures = mysort.util.dimReductionPCA(spikes, Rms.P.featureExtraction.nDims, [], 3*1000000);
                else
                    E{i,1}.spikeFeatures = Rms.S.spikeFeatures.X;
                end
            end

            R = struct();
            R.E = E;
        end
    end
end