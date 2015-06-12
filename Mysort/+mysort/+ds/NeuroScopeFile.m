classdef NeuroScopeFile < mysort.ds.PreFilteredDataSourceInterface & ...
                          mysort.ds.PreProcessedDataSourceInterface
    properties
        nC
        fname
        memmap
    end
    
    methods
        %------------------------------------------------------------------
        function self = NeuroScopeFile(fname, nC, s_per_sec, hpf, lpf)
            ME = mysort.ds.MultiElectrode([ones(nC,1) 17*(1:nC)'], ones(nC,1));
            [pathstr, name, ext] = fileparts(fname);
            prefilterfname  = fullfile(pathstr, [name '_prefilter_' num2str(hpf) '_' num2str(lpf) '.h5']);
            preprocessfname = fullfile(pathstr, [name '_preproc.h5']);
            filtererFactory_ = mysort.mea.FilterWrapper(hpf, lpf, s_per_sec, 10, 'butter');
            self = self@mysort.ds.PreFilteredDataSourceInterface(prefilterfname, '/sig', [], filtererFactory_, 1, 'NeuroScopeFile', s_per_sec, ME);
            self = self@mysort.ds.PreProcessedDataSourceInterface(preprocessfname, '/preproc', [], [], 'NeuroScopeFile', s_per_sec, ME);
            self.fname = fname;
            self.nC = nC;
            self.MultiElectrode.setDataSource(self);
            self.memmap = memmapfile(fname, 'Format', 'int16');
        end
        %------------------------------------------------------------------
        function self_ = copy(self)
            self_ = mysort.ds.NeuroScopeFile(self.fname, self.nC, self.samplesPerSecond, ...
                self.filterFactory.m_nHighPassFilter, self.filterFactory.m_nLowPassFilter);
        end
        %------------------------------------------------------------------
        function X = getRawData(self, timeIndex, channelIndex, varargin)
%             assert(max(timeIndex)<=size(self,1), 'Invalid timeIndex (too big)');
%             assert(min(timeIndex)>0, 'Invalid timeIndex (too small)');
%             assert(max(channelIndex)<=size(self,2), 'Invalid channel index (too big)');
%             assert(min(channelIndex)>0, 'Invalid channel (too small)');
            
            mi = (min(timeIndex) -1)*self.nC + 1;
            ma = (max(timeIndex) )*self.nC;
            X = self.memmap.data(mi:ma);
            X = reshape(X, self.nC, (ma-mi+1)/self.nC)';
            X = X(timeIndex-timeIndex(1)+1, channelIndex);
        end
        %------------------------------------------------------------------
        function X = getScaledData(self, timeIndex, channelIndex, varargin)
            P.progressDisplay = 'none'; % or 'console' or 'progressbar';
            P = mysort.util.parseInputs(P, '', varargin);   
            
            X = self.getRawData(timeIndex, channelIndex, ...
                'progressDisplay', P.progressDisplay);
         
%             X   = double(X)*self.lsb;
        end
        %------------------------------------------------------------------
        function X = getUnfilteredData(self, varargin)
            X = self.getScaledData(varargin{:});
        end
        %------------------------------------------------------------------
        function L = getNSamples_(self)
            L = length(self.memmap.data)/self.nC;
        end                   
    end
end