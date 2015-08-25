classdef CMOSMEA < mysort.ds.PreFilteredMultiSessionInterface & ...
                   mysort.ds.PreProcessedMultiSessionInterface
    properties
        P
        h5info
        fname
    end
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = CMOSMEA(fname, varargin)
            P.name = 'CMOSMEA';
            P.useFilter = 0;
            P.hpf = 500;
            P.lpf = 3000;
            P.filterOrder = 2;
            P.filterType = 'butter';
            P.filterName = '';
            P.prefilter = 0;
            P.preprocess = 0;
            P = mysort.util.parseInputs(P, varargin, 'error');
            
            assert(exist(fname, 'file')>0, sprintf('Could not find file: %s', fname));
            s_per_sec = hdf5read(fname, '/Sessions/Session0/sr');  
            %s_per_sec = s_per_sec(1);
            assert(length(s_per_sec)==1, 'samples per second is an array!?');
            
            h5info_ = hdf5info(fname);
            
            filterFactory_ = mysort.mea.FilterWrapper(P.hpf, ...
                P.lpf, s_per_sec, P.filterOrder, P.filterType);
            
            %myFolderName = fname(1:end-3);
            [pathstr, name] = fileparts(fname);
            myFolderName = fullfile(pathstr, name);
            if ~exist(myFolderName, 'file')
                mkdir(myFolderName);
            end
            filter_buffer_file_ = fullfile(myFolderName, ['filtered' ...
                num2str(P.filterOrder) '_' num2str(P.hpf)...
                '_' num2str(P.lpf) '.h5']);
            filter_buffer_file_h5info_ = [];
            if exist(filter_buffer_file_, 'file')
                filter_buffer_file_h5info_ = hdf5info(filter_buffer_file_);
            end
                
            preprocess_file_ = [filter_buffer_file_ 'preproc.h5'];
            preprocess_file_h5info_ = [];
            if exist(preprocess_file_, 'file')
                preprocess_file_h5info_ = hdf5info(preprocess_file_);
            end

            % init sessions
            nSessions = length(h5info_.GroupHierarchy(1).Groups(1).Groups);
            sessionList_ = mysort.mea.CMOSMEASession.empty(); 
            for i=1:nSessions
                sessionList_(i) = mysort.mea.CMOSMEASession(fname, h5info_,...
                    filter_buffer_file_, filter_buffer_file_h5info_, ...
                    preprocess_file_, preprocess_file_h5info_, ...
                    filterFactory_, i-1, P.useFilter);
            end
            [pathstr, name, ext] = fileparts(fname);
            sortingPath = fullfile(pathstr, name, 'sortings');
            self = self@mysort.ds.PreFilteredMultiSessionInterface(filter_buffer_file_, '/', filter_buffer_file_h5info_, filterFactory_, P.useFilter, P.name, s_per_sec, sessionList_);
            self = self@mysort.ds.PreProcessedMultiSessionInterface(preprocess_file_, '/', preprocess_file_h5info_, sortingPath, P.name, s_per_sec, sessionList_);
            for i=1:nSessions
                self.sessionList(i).parent = self;
            end
            self.fname = fname;
            self.h5info = h5info_;
            self.P = P;
            
            self.filterFactory = filterFactory_;
            
            if self.P.prefilter
               	self.prefilter();
            end

            if self.P.preprocess
                self.preprocess();
            end
        end

        %------------------------------------------------------------------
        function self2 = copy(self)
            self2 = mysort.mea.CMOSMEA(self.fname, 'hpf', self.P.hpf, 'lpf', self.P.lpf, 'filterOrder', self.P.filterOrder);
        end
        
        %------------------------------------------------------------------
        function X = getRawData(self, timeIndex, channelIndex, sessionIndex, varargin)
            if nargin < 4
                sessionIndex = self.activeSessionIdx;
                if nargin < 3
                    channelIndex = 1:self.size_(2);
                    if nargin < 2
                        timeIndex = 1:self.size_(1);
                    end
                end
            end    
            X = self.sessionList(sessionIndex).getRawData(timeIndex, channelIndex, varargin{:});
        end
        
        %
        function x = isBinaryFile(self)
            assert(~isempty(self.sessionList), 'Session list empty');
            s = self.sessionList(1);
            x = isa( s.h5matrix_raw, 'mysort.ds.binaryFileMatrix');
        end
        
        %------------------------------------------------------------------
        function W = getCutWaveforms(self, timepoints, cutleft, cutright, varargin)
            P.channels = [];
            P.session = self.activeSessionIdx;
            P.maxLoadSamples = 100000;
            P = mysort.util.parseInputs(P, varargin, 'error');

            if isempty(timepoints)
                W=[];
                return
            end
            
            L = cutright+cutleft+1;
            assert(L>0, 'cant cut negative epochs length!');
            nConnectedChannels = self.getNConnectedChannels();
            if isempty(P.channels);
                P.channels = 1:nConnectedChannels;
            else
                assert(max(P.channels)<=nConnectedChannels, 'channel index out of bounds')
            end
            
            if size(timepoints,1) == 1
                timepoints = timepoints';
            elseif size(timepoints,2) > 1
                error('timepoints must be a vector!')
            end
            
%             % group timepoints, so that all data of one group can be loaded
%             n = length(timepoints);
%             T = sortrows([timepoints; (1:n)']);
%             assert(any(sort(timepoints)-timepoints ~=0), 'timepoints must be sorted!');
            
            t1 = max(1, min(timepoints)-cutleft);
            t2 = min(self.size_(1, P.session), max(timepoints)+cutright);
            X = self.getData(t1:t2, P.channels, P.session)';
            epochs = [timepoints-cutleft timepoints+cutright]-t1 + 1;
%             eL = mysort.epoch.length(epochs);
%             cL = [0; cumsum(eL)];
%             IDX = zeros(1, cL(end));            
%             for i=1:size(epochs,1)
%                 IDX(cL(i)+1:cL(i+1)) = epochs(i,1):epochs(i,2);
%             end
%             Xs = X(:,IDX);
            W = mysort.epoch.extractWaveform(X, epochs);
        end


%         %------------------------------------------------------------------
%         function ME = getMultiElectrode(self, session_idx)
%             if nargin == 1
%                 session_idx = self.activeSessionIdx;
%             end
%             ME = self.sessionList(session_idx).getMultiElectrode();
%         end        
        %------------------------------------------------------------------
        function CL = getChannelList(self, session_idx)
            if nargin == 1
                session_idx = self.activeSessionIdx;
            end
            CL = self.sessionList(session_idx).getChannelList();
        end
        %------------------------------------------------------------------
        function CNr = getChannelNr(self, session_idx)
            if nargin == 1
                session_idx = self.activeSessionIdx;
            end
            CL = self.sessionList(session_idx).getChannelList();
            CNr = CL(CL(:,2)==1, 1);
        end        
        %------------------------------------------------------------------
        function N = getNConnectedChannels(self, session_idx)
            if nargin == 1
                session_idx = self.activeSessionIdx;
            end
            CL = self.sessionList(session_idx).getChannelList();
            N = sum(CL(:,2)==1);
        end        
        
        
        %------------------------------------------------------------------
        function sfnames = getSessionFilenames(self, varargin)
            sfnames = self.getSessionVar('sourceFname', varargin{:});
        end
        %------------------------------------------------------------------
        function messages = getSessionMessages(self, varargin)
            messages = self.getSessionVar('message', varargin{:});
        end 
    end
end