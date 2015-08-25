classdef Template < mysort.mea.Waveform
    properties
        nSourceSpikesPerElectrode
    end
    
    
    methods
        %------------------------------------------------------------------
        function self = Template(waveforms, cutLeft, samplesPerSecond, name, MultiElectrode, nSourceSpikesPerElectrode)
            % waveforms is time x channels
            if nargin == 1
                vars = {waveforms};
            else
                vars = {waveforms, cutLeft, samplesPerSecond, name, MultiElectrode};
            end
            self = self@mysort.mea.Waveform(vars{:});
            if nargin > 1
                assert(MultiElectrode.getNElectrodes() == length(nSourceSpikesPerElectrode), 'Electrodenumber does not match!');
                self.nSourceSpikesPerElectrode = nSourceSpikesPerElectrode;
            end
        end
        %------------------------------------------------------------------
        function fromStruct(self, S)
            fromStruct@mysort.mea.Waveform(self, S);
            self.nSourceSpikesPerElectrode = S.nSourceSpikesPerElectrode;
        end
        %------------------------------------------------------------------
        function S = toStruct(self)
            S = toStruct@mysort.mea.Waveform(self);
            S.nSourceSpikesPerElectrode = self.nSourceSpikesPerElectrode;
        end
        %------------------------------------------------------------------
        function cp = copy(self)
            S = self.toStruct();
            cp = mysort.mea.Template(S);
        end
        %------------------------------------------------------------------
        function T = getTemplateWithHighestChannels(self, nMax)
            if size(self.waveforms,2) <= nMax
                T = self;
                return
            end
            maxChans = max(abs(self.waveforms), [], 1)';
            maxChans(isnan(maxChans)) = -1;
            maxChans = sortrows([maxChans (1:size(self.waveforms,2))']);

            plotChanIdx = sort(maxChans(end-nMax+1:end,2));            
            T = self.getTemplate4ElectrodeIdx(plotChanIdx);
        end
        %------------------------------------------------------------------
        function T = getTemplate4ElectrodeIdx(self, idx)
            T = mysort.mea.Template(self.waveforms(:,idx),...
                self.cutLeft, self.samplesPerSecond, self.name, ...
                self.MultiElectrode.getSubElectrode4ElIdx(idx), self.nSourceSpikesPerElectrode(idx));
        end  
        %------------------------------------------------------------------
        function n = getNSourceWaveforms4ElectrodeIdx(self, idx)
            n = self.nSourceSpikesPerElectrode(idx);
        end
        %------------------------------------------------------------------
        function n = getNSourceWaveforms4MultiElectrode(self, ME)
            requestedElectrodes = ME.electrodeNumbers;
            n = zeros(1, length(requestedElectrodes));
            [C ia ib] = intersect(self.MultiElectrode.electrodeNumbers, requestedElectrodes);
            n(ib) = self.getNSourceWaveforms4ElectrodeIdx(ia);
        end
        
        
        %------------------------------------------------------------------
        function cp = merge(self, T)
            assert(self.samplesPerSecond == T.samplesPerSecond, 'Sampling rates must be identical!');
            assert(self.cutLeft == T.cutLeft, 'cutLeft has to be identical');
            if self.MultiElectrode.getNElectrodes() > 0
                assert(self.MultiElectrode.hasElectrodeNumbers(), 'MultiElectrode must have electrode numbers to be mergable!');
            end
            if T.MultiElectrode.getNElectrodes() > 0
                assert(T.MultiElectrode.hasElectrodeNumbers(), 'MultiElectrode must have electrode numbers to be mergable!');
            end
            cp = self.copy();
            [cp.MultiElectrode mergedIdx] = self.MultiElectrode.merge(T.MultiElectrode);
            
            for el = 1:T.MultiElectrode.getNElectrodes()
                idx = mergedIdx(el);
                if size(self.waveforms,2) >= idx
                    % the waveform of that channel already exists
                    n1 = cp.nSourceSpikesPerElectrode(idx);
                    n2 = T.nSourceSpikesPerElectrode(el);
                    cp.waveforms(:,idx) = (n1*cp.waveforms(:,idx) + n2*T.waveforms(:,el))/(n1+n2);
                    cp.nSourceSpikesPerElectrode(idx) = n1+n2;
                else
                    % create new channel
                    cp.waveforms(:,idx) = T.waveforms(:,el);
                    cp.nSourceSpikesPerElectrode(idx) = T.nSourceSpikesPerElectrode(el);
                end
            end
        end
    end
end