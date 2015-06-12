classdef CMOSMEASession < mysort.ds.PreFilteredDataSourceInterface & ...
                          mysort.ds.PreProcessedDataSourceInterface
    properties
        parent
        fname
        sourceFname
        message
        h5matrix_raw
        session_idx
        session_str
        size_buffer
        CL
        connected_channels
        lsb
    end
    
    methods
        %------------------------------------------------------------------
        function self = CMOSMEASession(fname, h5info_var, prefilterfname, prefh5info, preprocessfname, preproch5info, filterFactory, session_idx, useFilter)
            session_str_ = sprintf('/Sessions/Session%d/', session_idx);
            s_per_sec = hdf5read(fname, [session_str_ 'sr']);   
            % HACK !!! (this is ok, was a bug in the old ntk converter from
            % Jan.
            s_per_sec = s_per_sec(1);
            
            % build the multielectrode
            % CHECK IF THIS IS A NEW FILE VERSION WITH THE CHANNEL LIST
            % SPLIT UP
            try 
                % THIS IS MUCH FASTER
                connectedChannels = hdf5read(fname, [session_str_ 'channel_connected']);
                nC_ = length(connectedChannels);
                connectedChannels = connectedChannels == 1;
                assert(any(connectedChannels), 'None of the channels was connected!')
                cpx = double(hdf5read(fname, [session_str_ 'channel_posx']));
                cpy = double(hdf5read(fname, [session_str_ 'channel_posy']));
                cnr = double(hdf5read(fname, [session_str_ 'channel_nr']));
                cpx = cpx(connectedChannels==1);
                cpx = cpx(:);
                cpy = cpy(connectedChannels==1);
                cpy = cpy(:);
                ME = mysort.ds.MultiElectrode([cpx cpy], cnr);
            catch
                % Seems to be the old format
                % WHY EVER, BUT THIS READ CALL IS SLOW LIKE HELL
                CL__ = hdf5read(fname, [session_str_ 'channel_list']);
                CL__ = get(CL__, 'Data');
                nC = size(CL__,1);
                nV = length(CL__{1});
                CL_ = zeros(nC, nV);
                for i=1:nC
                    CL_(i,:) = cellfun(@double, CL__{i});
                end
                connectedChannels = CL_(:,2)==1;
                CL_(~connectedChannels, :) = [];   % remove unconnected channels
                CL_(:,3:4) = CL_(:,3:4)/1000; % convert to micro meter
                ME = mysort.ds.MultiElectrode(CL_(:,3:4), CL_(:,5));
            end
            
            
            self = self@mysort.ds.PreFilteredDataSourceInterface(prefilterfname, [session_str_ 'sig'], prefh5info, filterFactory, useFilter, 'CMOSMEASession', s_per_sec, ME);
            self = self@mysort.ds.PreProcessedDataSourceInterface(preprocessfname, [session_str_ 'preproc'], preproch5info, [], 'CMOSMEASession', s_per_sec, ME);
            self.fname = fname;
            self.session_idx = session_idx;
            self.session_str = session_str_;
            self.size_buffer = [];
            
            %find out whether sig is data of a string if a binary file name:
            fileInfo = h5info(self.fname, [self.session_str 'sig']);
            is_int = strcmp(fileInfo.Datatype.Class, 'H5T_INTEGER');
            is_str = strcmp(fileInfo.Datatype.Class, 'H5T_STRING') ;
            assert(is_int || is_str, 'sig has to be in the format of an integer or a string');
            
            if is_str
                binFileName = cell2mat( h5read(self.fname, [self.session_str 'sig']) );
                
                [pathstr,name,ext] = fileparts(self.fname);
                binFile = fullfile(pathstr, binFileName);
                
                binDims = h5read(self.fname, [self.session_str 'bin_dims']);
                assert( exist(binFile, 'file') == 2, ['Task aborted: binary file ' binFile ' not found!']);
                % todo: naming not good!
                self.h5matrix_raw =  mysort.ds.binaryFileMatrix(binFile, binDims);
            else
                self.h5matrix_raw = mysort.h5.matrix(self.fname, [self.session_str 'sig'], true);
            end
            
            try
                self.sourceFname = get(hdf5read(self.fname, [self.session_str 'fileName']), 'data');
            catch
                self.sourceFname = [];
            end
            self.message = 'unkown';
            if mysort.h5.exist(self.fname, [self.session_str 'message'], h5info_var)
                self.message = get(hdf5read(self.fname, [self.session_str 'message']), 'data');
            end
            self.MultiElectrode.setDataSource(self);

            % get lsb:
            self.lsb = self.getLSB();
            self.connected_channels = find(connectedChannels);
        end

        %------------------------------------------------------------------
        function gain = getGain(self)
            gain = hdf5read(self.fname, [self.session_str 'gain']);
            assert(length(gain) == 4, 'Gain must have 4 values!'); 
            try
                gainmultiplier = hdf5read(self.fname, [self.session_str 'filter/gainmultiplier']);
            catch
                gainmultiplier = 1;
            end
            gain(1) = gain(1)*single(gainmultiplier);
        end
        %------------------------------------------------------------------
        function lsb = getLSB(self)
            try
                ntk.adc_resolution = hdf5read(self.fname, [self.session_str 'adc_resolution']);
            catch
                ntk.adc_resolution = 8;
            end
            try
                ntk.adc_range = hdf5read(self.fname, [self.session_str 'adc_range']);
            catch
                ntk.adc_range = 2.9;
            end            
            gain = self.getGain();
            if ~isempty(gain) && any(gain~=1)
                lsb = ntk.adc_range/(2^ntk.adc_resolution-1)/gain(1)*1000000;
            else
                lsb = 1;
            end
        end
        %------------------------------------------------------------------
        function CL = getChannelList(self)
            if isempty(self.CL)
                CL_ = hdf5read(self.fname, [self.session_str 'channel_list']);
                CL_ = get(CL_, 'Data');
                nC = size(CL_,1);
                nV = length(CL_{1});
                self.CL = zeros(nC, nV);
                for i=1:nC
                    self.CL(i,:) = cellfun(@double, CL_{i});
                end
            end
            CL = self.CL;
        end        
        
        %------------------------------------------------------------------
        function wf = getWaveform_(self, varargin)
            wf = self.h5matrix_raw.getWaveform_(varargin{:});
        end

        
        %------------------------------------------------------------------
        function X = getRawData(self, timeIndex, channelIndex, varargin)
            P.channels = []; %set channels to 'all' to also get unconnected channels
            P = mysort.util.parseInputs(P, '', varargin);
            
            if isempty(P.channels) || strcmp(P.channels, 'connected')
                assert(max(channelIndex)<=length(self.connected_channels), 'channel index out of bounds')
                channelIndex = self.connected_channels(channelIndex);
            end
            % channel idx is first dimension !!!
            X = self.h5matrix_raw(timeIndex, channelIndex);
        end
        %------------------------------------------------------------------
        function X = getScaledData(self, timeIndex, channelIndex, varargin)
%             P.progressDisplay = 'none'; % or 'console' or 'progressbar';
%             P = mysort.util.parseInputs(P, '', varargin);   
            
            channelIndex = self.connected_channels(channelIndex);
            tmp = self.h5matrix_raw(timeIndex, channelIndex);
            X = double(tmp)*self.lsb;
            % from ntk_stream_process (line 130):
            % y.gain=ntk.gain1*ntk.gain2*ntk.gain3;
            % y.lsb=ntk.adc_range/(2^ntk.adc_resolution-1)/y.gain*1000000;
            % y.range=y.lsb*(2^ntk.adc_resolution-1);
            
            % from ntk_stream_process (line 291):
            % y.sig=double(ntk.data(keep_channels,:)')*y.lsb;
             
            % get gain
            %range = lsb*(2^ntk.adc_resolution-1);
        end
        %------------------------------------------------------------------
        function X = getUnfilteredData(self, varargin)
            X = self.getScaledData(varargin{:});
        end
        %------------------------------------------------------------------
        function L = getNSamples_(self)
%             if isempty(self.size_buffer)
%                 CL = self.getChannelList();
%                 cols = sum(CL(:,2));
%                 path_str = [self.session_str 'sig'];
%                 [bFound session_data Dims] = mysort.h5.exist(self.fname, path_str);
%                 rows = Dims(1);
%                 channel = Dims(2);
%                 assert(bFound, 'The group hierarchy inside this H5 file is not consistent with the current implementation. Dimensions of Data (sig) cannot be established!');
%                 assert(size(CL,1) == channel, 'The number of rows in the H5 file is not the same as the number of channels in the channel list!');
%                 self.size_buffer = [rows cols];
%             else
%                 rows = self.size_buffer(1);
%                 cols = self.size_buffer(2);
%             end
            L = size(self.h5matrix_raw,1);
        end   
        %------------------------------------------------------------------
        function Cest = getCovest(self, varargin)
            if self.bCovIsPrecalculated
                Cest = self.getCovestFromBufferFile(varargin{:});
            else
                Cest = [];
                % Try to get Cest from some session before us
                if self.session_idx > 0 && ~isempty(self.parent)
                    ME1 = self.MultiElectrode;
                    i = 0;
                    while i<self.session_idx && isempty(Cest)
                        ME2 = self.parent.sessionList(i+1).MultiElectrode;
                        if length(ME1.electrodeNumbers) == length(ME2.electrodeNumbers) && ...
                           ~any(ME1.electrodeNumbers ~= ME2.electrodeNumbers)
                            Cest = self.parent.sessionList(i+1).getCovest(varargin{:});
                        end
                        i = i+1;
                    end
                end
                % If Cest could not be retrieved from someone else,
                % calculate with super method
                if isempty(Cest)
                    Cest = getCovest@mysort.ds.ExtendedDataSourceInterface(self, varargin{:});
                end
                CestS = Cest.toStruct();
                save(self.preprocCov, 'CestS');
                self.bCovIsPrecalculated = true;                
            end
        end                  
    end
end