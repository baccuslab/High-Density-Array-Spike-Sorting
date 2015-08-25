classdef WaveformMatrix < mysort.ds.WaveformDataSourceInterface
    properties
        ME
        X 
        cutLeft
        cutLength
    end
      
    methods
        %------------------------------------------------------------------
        function self = WaveformMatrix(X, varargin)
            P.name = 'WaveformMatrix';
            P.samplesPerSecond = [];
            P.nC = 1;
            P.cutLeft = 0;
            P.multiElectrode = [];
            P = mysort.util.parseInputs(P, varargin);
            if isempty(P.multiElectrode)
                P.multiElectrode = mysort.ds.MultiElectrode(1:P.nC);
            end
            self = self@mysort.ds.WaveformDataSourceInterface(P.name, P.samplesPerSecond, P.multiElectrode);

            assert(isnumeric(X), 'data must be a matrix!');
            assert(ndims(X)==2, 'data must be a matrix!'); 
            %cutLength_ = size(X,2)/self.ME.getNElectrodes();
            cutLength_ = size(X,2)/self.MultiElectrode.getNElectrodes();
            assert(cutLength_ == round(cutLength_), 'Number of channels does not match data matrix!')
            self.X = X;
            self.cutLeft = P.cutLeft;
            self.cutLength = cutLength_;
        end
        %------------------------------------------------------------------
        function wfs = getWaveform(self, rowIdx, cutLeft, cutLength, channelidx)
            if nargin < 3 || isempty(cutLeft)
                cutLeft = self.cutLeft;
            end
            if nargin < 4 || isempty(cutLength)
                cutLength = self.cutLength;
            end
            assert(cutLength>0 && cutLength<=self.cutLength, 'cutLength out of bounds!')
            assert(cutLeft <= self.cutLeft, 'cutLeft out of bounds!')
            assert(cutLength-cutLeft <= self.cutLength, 'cutLeft + cutLength combination is out of bounds!');
            wfs = self.X(rowIdx,:);
            nC = self.ME.getNElectrodes();
            if nargin == 5 && ~isempty(channelidx)
                % we need to subselect channels
                wfs = mysort.wf.vSubChanSel(self.X(rowIdx,:), self.ME.getNElectrodes(), channelidx); 
                nC = length(channelidx);
            end
            if ~(cutLength == self.cutLength && cutLeft == self.cutLeft)
                % we need to cut subset
                t1 = self.cutLeft - cutLeft +1;
                t2 = t1+cutLength-1;
                tidx = t1:t2;
                wfs = mysort.wf.vSubsel(wfs, nC, tidx);
            end
        end
    end
end