
classdef ConcatenatedWaveform < mysort.plot.SingleAxesPlot
    properties (Constant)
    end
    
    properties
    end
    
    methods (Abstract)
    end
    
    methods 
        %%% ------------------------------------------------------ 
        function self = ConcatenatedWaveform(waveforms, varargin)
            self = self@mysort.plot.SingleAxesPlot(varargin{:});
            self.P.nC = 1;
            self.P.Id = 1;
            self.P.linewidth = 2;
            self.P.colorWithId = false;
            self.P.NoXTickLabel = 0;
            self.P.channelIDs = []; 
            self.P.plotMean = 1;
            self.P = mysort.util.parseInputs(self.P, 'ConcatenatedWaveform', varargin, 1);

            nC = self.P.nC;
            if ~isempty(self.P.channelIDs)
                if nC == 1
                    nC = length(self.P.channelIDs);
                end
                assert(nC==length(self.P.channelIDs), 'You must provide as many channel IDs as there are channels!');
            end
            nS = size(waveforms,1);
            Tf = size(waveforms,2)/nC;
            assert(round(Tf) == Tf, 'nC is not correct!');
            
            xRange = [];
            plotWavs = nan(nS, Tf*nC + nC-1);
            for c=1:nC
                xRange = [xRange (c-1)*Tf + (1:Tf) nan];
                plotWavs(:,(c-1)*(Tf+1)+1:c*(Tf+1)-1) = [waveforms(:, (c-1)*Tf+1:c*Tf)];
            end
            xRange = xRange(1:end-1);
            if self.P.colorWithId
                plot(repmat(xRange, nS, 1)', plotWavs', 'linewidth', self.P.linewidth, 'color', mysort.plot.vectorColor(self.P.Id));
            else
                plot(repmat(xRange, nS, 1)', plotWavs', 'linewidth', self.P.linewidth);
            end
            hold on
            if self.P.plotMean
                plot(repmat(xRange, nS, 1)', mean(plotWavs), 'linewidth', 2, 'color', 'k');
            end
            axis tight
            set(self.P.ax,'xlim', [1 Tf*nC + 10]);
            set(self.P.ax,'TickLength', [0.0001; 0.02]);
            if self.P.NoXTickLabel
                set(self.P.ax, 'XTickLabel', []);
            end
            plot(Tf*nC+5, 0, '.', 'markerSize', 25, 'color', mysort.plot.vectorColor(self.P.Id));
            box off
            mysort.plot.verticalLines(((1:nC-1)*Tf)+.5, [], 'k:');
            ylabel(sprintf('%d | %d', self.P.Id, nS));
            if ~isempty(self.P.channelIDs)
                pos = get(self.P.ax, 'position');
                fh = get(self.P.ax, 'parent');
                for i=1:nC
                    p = pos;
                    p(1) = pos(1) +  (i-1)/nC*pos(3);
                    p(2) = pos(2) + .99*pos(4);
                    p(3:4) = [.1 .04];
                    if isnumeric(self.P.channelIDs(i))
                        str = ['channel: ' num2str(self.P.channelIDs(i))];
                    else
                        str = ['channel: ' self.P.channelIDs(i)];
                    end
                    annotation(fh, 'textbox', p, 'String', str, 'LineStyle', 'none');
                end
            end
        end
    end
end