classdef ntkFileWrapper < handle
    properties (SetAccess=private)
        
    end
    properties
        P
        filename
        ntk
        bNtkIsInitialized
        MultiElectrode
        MultiElectrodeFull
        isConnectedElectrode
    end
    
    methods
        %%% ----------------CONSTRUCTOR---------------------------
        function self = ntkFileWrapper(fname, varargin)
            self.P.chunkSize = 2e6;
            self.P = mysort.util.parseInputs(self.P, varargin, 'error');
            self.filename = fname;
            
            self.init();
        end
        %%% ------------------------------------------------------
        function self = init(self, varargin)
            siz = 10;
            self.ntk = initialize_ntkstruct(self.filename, 'nofilters');
            [ntk_raw self.ntk] = ntk_load(self.ntk, siz);
            nC     = size(ntk_raw.sig, 2);
            electrodePositions = [ntk_raw.x' ntk_raw.y']; 
            electrodeNumbers = ntk_raw.channel_nr';
            self.MultiElectrode = mysort.ds.MultiElectrode(electrodePositions, electrodeNumbers);
            self.bNtkIsInitialized = false;
%             self.isConnectedElectrode = zeros(1, size(electrodePositionsAll,1));
%             for i=1:nCFull
%                 self.isConnectedElectrode(i) = ~isempty(self.ntk.channels{i,1}.els);
%             end
%             self.MultiElectrode = mysort.ds.MultiElectrode(electrodePositionsAll(self.isConnectedElectrode==1),...
%                 electrodeNumbersAll(self.isConnectedElectrode==1));
        end
        %%% ------------------------------------------------------
        function X = getAllFilteredData(self)
            X = self.prefilter(double(self.getAllRawData()));
        end
        %%% ------------------------------------------------------
        function X = getNextChunkFilteredData(self)
            X = self.prefilter(double(self.getNextChunkRawData()));
        end        
        %%% ------------------------------------------------------
        function X = prefilter(self, X)
            tic
            X = bsxfun(@minus, X, mean(X,1));
            toc
            tic
            X = conv2(X, mysort.mea.filter_gpu()','same');
            toc
        end
        %%% ------------------------------------------------------
        function X = getNextChunkRawData(self)
            siz = self.P.chunkSize;
            if ~self.bNtkIsInitialized
                self.ntk = initialize_ntkstruct(self.filename, 'nofilters');
                self.bNtkIsInitialized = true;
            end
            if self.eof()
                warning('End of file reached!')
                X = [];
                return
            end
            [high_density_data self.ntk] = ntk_load(self.ntk, siz);
            X = double(high_density_data.sig);
        end         
        %%% ------------------------------------------------------
        function b = eof(self)
            if ~self.bNtkIsInitialized
                b = false;
                return
            end
            b = self.ntk.eof;
        end
        %%% ------------------------------------------------------
        function X = getAllRawData(self)
            siz = 1e9;
            ntk = initialize_ntkstruct(self.filename, 'nofilters');
            [high_density_data ntk] = ntk_load(ntk, siz);
            X = double(high_density_data.sig);
            % LOAD DATA
            while ~ntk.eof
                [high_density_data ntk] = ntk_load(ntk, siz);
                X = [X; double(high_density_data.sig)];
            end       
        end   
    end
end