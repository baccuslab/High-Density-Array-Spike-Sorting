classdef Website < mysort.ds.PreFilteredDataSourceInterface &...
                   mysort.ds.ExtendedDataSourceInterface
    properties
        fname
        h5matrix_raw
    end
    
    methods
        %------------------------------------------------------------------
        function self = Website(fname, hpf, lpf)
            samplesPerSecond   = double(hdf5read(fname, '/samples_per_second'));
            electrodePositions = double(hdf5read(fname, '/electrode_pos'))';
            electrodeNumbers = 1:size(electrodePositions,1);
            ME = mysort.ds.MultiElectrode(electrodePositions, electrodeNumbers); 
            prefilterBufferFilename = [fname '_pref'];
            prefilterBufferH5path = '/data';
            h5info = [];
            FF = mysort.mea.FilterWrapper(hpf, lpf, samplesPerSecond, 6, 'butter');
            self = self@mysort.ds.PreFilteredDataSourceInterface(prefilterBufferFilename, prefilterBufferH5path, h5info, FF, 'BenchmarkDS', samplesPerSecond, ME);
            self = self@mysort.ds.ExtendedDataSourceInterface('BenchmarkDS', samplesPerSecond, ME);
            self.h5matrix_raw = mysort.h5.matrix(fname, '/data');
            self.fname = fname;
        end
        %------------------------------------------------------------------
        function X = getUnfilteredData(self, timeIndex, channelIndex, varargin)
            % channel idx is first dimension !!!
%             X = self.h5matrix_raw(timeIndex, channelIndex);
            X = double(self.h5matrix_raw(channelIndex, timeIndex))';
        end

        %------------------------------------------------------------------
        function L = getNSamples_(self)
%             L = size(self.h5matrix_raw,1);
            L = size(self.h5matrix_raw,2);
        end   
    end
end