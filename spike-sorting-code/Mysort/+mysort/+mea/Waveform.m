classdef Waveform < handle
    properties
        name
        waveforms % tensor, "Length x Electrodes x nWaveforms"
        samplesPerSecond
        MultiElectrode
        cutLeft
    end
    
    
    methods
        %------------------------------------------------------------------
        function self = Waveform(waveforms, cutLeft, samplesPerSecond, name, MultiElectrode)
            if ischar(waveforms) && ischar(samplesPerSecond)
                fname = waveforms;
                h5path = samplesPerSecond;
                self.loadFromFile(fname, h5path)
                return                
            end
            if isstruct(waveforms)
                self.fromStruct(waveforms);
                return
            end
            self.waveforms = waveforms;
            self.cutLeft = cutLeft;
            self.samplesPerSecond = samplesPerSecond;
            self.name = name;
            assert(MultiElectrode.getNElectrodes() == size(waveforms,2), 'Waveforms must be in tensor form: time x channel x idividuals!');
            self.MultiElectrode = MultiElectrode;
        end
        %------------------------------------------------------------------
        function fromStruct(self, S)
            self.name = S.name;
            self.samplesPerSecond = S.samplesPerSecond;
            self.cutLeft = S.cutLeft;
            self.MultiElectrode = mysort.ds.MultiElectrode(S.MultiElectrode);
            nC = self.MultiElectrode.getNElectrodes();
            self.waveforms = mysort.wf.v2t(S.waveforms, nC);
        end
        %------------------------------------------------------------------
        function S = toStruct(self)
            S.name = self.name;
            S.samplesPerSecond = self.samplesPerSecond;
            S.waveforms = mysort.wf.t2v(self.waveforms);
            S.cutLeft = self.cutLeft;
            S.MultiElectrode = self.MultiElectrode.toStruct();
            S.date_ = date();
            S.version = 1;
            S.readme = 'This is a struct derived from the class mea.Waveform. Dont edit if you dont know what you are doing.';
        end        
        %------------------------------------------------------------------
        function save2File(self, fname, h5path)
            S = self.toStruct();
            mysort.h5.recursiveSave(fname, S, h5path);
        end
        %------------------------------------------------------------------
        function loadFromFile(self, fname, h5path)
            S = mysort.h5.recursiveLoad(fname, h5path);
            self.fromStruct(S);
        end        
        
        %------------------------------------------------------------------
        function wf = getWaveform4MultiElectrode(self, ME)
            requestedElectrodes = ME.electrodeNumbers;
            nC =  length(requestedElectrodes);
            Tf = size(self.waveforms,1);
            wf = zeros(nC, Tf);
            for i=1:length(requestedElectrodes)
                el = requestedElectrodes(i);
                idx = find(self.MultiElectrode.electrodeNumbers == el,1);
                if ~isempty(idx)
                    wf(i,:) = self.waveforms(:,idx);
                end
            end
        end
    end
end