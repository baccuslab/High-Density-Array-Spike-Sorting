classdef preprocessGroundTruth < ana.botmpaper.harrisAnalysis.analysis
    properties
        PREC
    end
    
    methods
        %------------------------------------------------------------------
        function self = preprocessGroundTruth(PREC)
            % Input argument PREC is optional but can be an
            % ana.botmpaper.harrisAnalysis.preprocess object
            self = self@ana.botmpaper.harrisAnalysis.analysis('PreprocessGT');
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
            
            activeIdxI = find(R.IntraAlignmentMask);
            nC = size(R.tExtraSpikesAlignedIntra,2);
            TfE = size(R.tAlignedExtraSpikesUpsampled,1);
            activeIdxE = [];
            for i=1:nC
                activeIdxE = [activeIdxE find(R.ExtraAlignmentMask(i,:))+(i-1)*TfE];
            end         
            
            % Produce Figure showing intracellular and extracellular spikes
            % aligned on INTRACELLULAR trace
            fh(p+1) = mysort.plot.figure('w', 1200, 'h', 1000); p=p+1;
                ah(1) = subplot(4,1,1);
                plot(1:size(R.tAlignedIntraSpikesUpsampled,1), squeeze(R.tAlignedIntraSpikesUpsampled))
                hold on
                plot(activeIdxI, zeros(1, length(activeIdxI)), 'rx');
                title('Upsamples and aligned waveforms with alignment mask');
                
                ah(2) = subplot(4,1,2);
                plot((1:size(R.tIntraSpikesAlignedIntra,1)), squeeze(R.tIntraSpikesAlignedIntra));
                title('Recut intracellular waveforms');            
                
                ah(3) = subplot(4,1,3:4);
                mysort.plot.waveformsVertical(R.tExtraSpikesAlignedIntra, 'IDs', ones(size(R.tExtraSpikesAlignedIntra,3),1),...
                    'axesHandle', ah(3), 'plotMedian', 1)
                title('Recut extracellular waveforms'); 
                linkaxes(ah(2:3), 'x');
                mysort.plot.figureTitle([name '- Alignment: Intra']);
            figureNames{p} = 'IntracellularAlignment';
            
            % Produce Figure showing intracellular and extracellular spikes
            % aligned on EXTRACELLULAR trace
            fh(p+1) = mysort.plot.figure('w', 1200, 'h', 1000); p=p+1;
                ah(1) = subplot(4,1,1);
                vS = mysort.wf.t2v(R.tAlignedExtraSpikesUpsampled(:,:,1:8:end));
                mysort.plot.spikes(vS, 'nC', nC, 'ax', ah(1))
                hold on
                plot(activeIdxE, zeros(1, length(activeIdxE)), 'rx');
                title('Upsamples and aligned waveforms with alignment mask');
                
                ah(2) = subplot(4,1,2);
                plot((0:size(R.tIntraSpikesAlignedExtra,1)-1)-R.gtAlignmentCutLeft, squeeze(R.tIntraSpikesAlignedExtra));
                title('Recut intracellular waveforms');            
                
                ah(3) = subplot(4,1,3:4);
                mysort.plot.waveformsVertical(R.tExtraSpikesAlignedExtra, 'IDs', ones(size(R.tExtraSpikesAlignedExtra,3),1),...
                    'axesHandle', ah(3), 'xOffset', -R.gtAlignmentCutLeft-1, 'plotMedian', 1)
                title('Recut extracellular waveforms'); 
                linkaxes(ah(2:3), 'x');
                mysort.plot.figureTitle([name '- Alignment: Extra']);
            figureNames{p} = 'ExtracellularAlignment';   
                
        end
    end
    methods(Access=protected)
        %------------------------------------------------------------------
        function R = runTrial_(self, name)
            R = self.PREC.getAllVariables(name);
            cutLeft = 30;
            if isempty(R) || isempty(R.tgdf)
                R = 'No ground truth available';
                return
            end
            trueintraspikes = mysort.epoch.extractWaveform(R.I, [R.tgdf(:,2)-cutLeft R.tgdf(:,2)+150]);
            nC = size(R.Xfil,1);
            % remove baseline on intracellular spikes
            baseline = mean(mean(trueintraspikes(:,[1 end])));
            trueintraspikes = trueintraspikes-baseline;
            
            tI = mysort.wf.v2t(trueintraspikes, 1);
            resampleFactor = 5;
            [tIA iTau_up meanMaskI] = mysort.wf.tAlignOnUpsampleMean(tI, 'upsample', resampleFactor,...
                'maxIter', 100, ...
                'downsample', 0, ...
                'useMedian', 0, ...
                'initAlignment', [], ... % DO NOT FORGET THIS !!!
                'restrictToNMaximalValues', 40,...
                'meanMaskMinLen', 1.5*resampleFactor,...
                'maxShiftPerIter', 15);
            meanAI = mean(squeeze(tIA),2);
            [m maxIdxIA] = max(meanAI);
            
            
            % make the index of the intracellular maximum the time of the
            % spike
            tgdfI = R.tgdf;
            tgdfI(:,2) = tgdfI(:,2) - iTau_up/resampleFactor + maxIdxIA/resampleFactor - cutLeft;
            
            trueextraspikesAI = mysort.epoch.extractWaveform(R.Xfil, round([tgdfI(:,2)-cutLeft tgdfI(:,2)+40]));
            trueintraspikesAI = mysort.epoch.extractWaveform(R.I, round([tgdfI(:,2)-cutLeft tgdfI(:,2)+40]));
            tEAI = mysort.wf.v2t(trueextraspikesAI, nC);
            tIAI = mysort.wf.v2t(trueintraspikesAI, 1);
            
            [tEA eTau_up meanMaskE] = mysort.wf.tAlignOnUpsampleMean(tEAI, 'upsample', resampleFactor,...
                'maxIter', 100, ...
                'downsample', 0, ...
                'useMedian', 1, ...
                'initAlignment', [], ... % DO NOT FORGET THIS !!!
                'restrictToNMaximalValues', 40,...
                'meanMaskMinLen', 1.5*resampleFactor,...
                'maxShiftPerIter', 1);
            meanAE = mean(squeeze(mean(tEA,2)),2);
            [m maxIdxEA] = min(meanAE);


            % make the index of the extracellular minimum the time of the
            % spike
            tgdfE = tgdfI;
            tgdfE(:,2) = tgdfE(:,2) - eTau_up/resampleFactor + maxIdxEA/resampleFactor -cutLeft ;
            
            cutLeft = 10;
            trueSpikesAE = mysort.epoch.extractWaveform([R.I; R.Xfil], round([tgdfE(:,2)-cutLeft tgdfE(:,2)+40]));
            tSAE = mysort.wf.v2t(trueSpikesAE, nC+1);
            tEAE = tSAE(:,2:end,:);
            tIAE = tSAE(:,1,:);

            R = struct();
            R.tgdfI = tgdfI;
            R.tgdfE = tgdfE;
            R.tExtraSpikesAlignedIntra = tEAI;
            R.tExtraSpikesAlignedExtra = tEAE;
            R.tIntraSpikesAlignedIntra = tIAI;
            R.tIntraSpikesAlignedExtra = tIAE;
            R.gtAlignmentCutLeft = cutLeft;
            R.alignmentResampleFactor = resampleFactor;

            R.tAlignedIntraSpikesUpsampled = tIA;
            R.tAlignedExtraSpikesUpsampled = tEA;
            R.IntraAlignmentMask = meanMaskI;
            R.ExtraAlignmentMask = meanMaskE;  
        end
    end
end