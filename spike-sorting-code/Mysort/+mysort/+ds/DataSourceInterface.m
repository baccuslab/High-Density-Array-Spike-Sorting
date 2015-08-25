classdef DataSourceInterface < mysort.ds.WaveformDataSourceInterface
    properties
        activeChannels
        fullMultiElectrode
    end
    
    methods(Abstract)
        getData_(self, idx1, idx2)
        getNSamples_(self)
    end
    
    methods
        %------------------------------------------------------------------
        function self = DataSourceInterface(varargin)
            self = self@mysort.ds.WaveformDataSourceInterface(varargin{:});
        end
        
        %------------------------------------------------------------------
        function varargout = subsref(self,S)
            bIsObjectArray = length(self)~=1;
            bIsIndexAccess = length(S)==1 && strcmp(S.type, '()');
            %             bIsMemberAccess = length(S) > 1 && strcmp(S(1).type, '.');
            %             bIsMethodCall   = bIsMemberAccess && ismethod(self, S(1).subs);
            %             bIsFieldAccess = bIsMemberAccess && ~bIsMethodCall && ismember(S(1).subs, fieldnames(self));
            %
            %             if bIsFieldAccess
            %                 varargout =
            %             end
            %            bIsFunctionAccessFromArrayAccess = (length(S)>1 && strcmp(S(1).type, '()') && strcmp(S(2).type, '.'));
            if ~bIsObjectArray && bIsIndexAccess && length(S.subs)==1
                assert(S.subs{1}==1, 'indexing out of bounds!')
                varargout{1} = self;
                return
            end
            if bIsObjectArray || ~bIsIndexAccess
                % for the '.' subsref call the standard one
                [varargout{1:nargout}] = builtin('subsref', self, S);
                return
            end
            assert(strcmp(S.type, '()'), 'Only () is implemented!');
            assert(length(S.subs) > 1, '(x) is not implemented!')
            varargout{1} = self.getData(S.subs{:});
        end
        %------------------------------------------------------------------
        function varargout = size(self,varargin)
            dims = self.getDims();
            varargout = matlabfilecentral.parseSize.parseSize(dims,nargout,varargin{:});
        end
        %------------------------------------------------------------------
        function out = end(self, k, n)
            out = self.size(k);
        end
        %------------------------------------------------------------------
        function dims = getDims(self)
            dims = [self.getNSamples() self.getNChannels()];
        end
        
        %------------------------------------------------------------------
        function L = getNSamples(self)
            L = self.getNSamples_();
        end
        %------------------------------------------------------------------
        function n = getNChannels(self)
            n = self.MultiElectrode.getNElectrodes();
        end
        
        %------------------------------------------------------------------
        function self = transpose(self)
            error('not implemented')
        end
        %------------------------------------------------------------------
        function self = ctranspose(self)
            error('not implemented')
        end
        %------------------------------------------------------------------
        function restrictToChannels(self, channelidx)
            if isempty(self.fullMultiElectrode)
                % Make a copy of full ME before setting the sub ME as
                % active
                self.fullMultiElectrode = self.MultiElectrode;
            end
            if nargin == 1 || isempty(channelidx)
                % reset the full ME to be active
                self.MultiElectrode = self.fullMultiElectrode;
                self.activeChannels = [];
            else
                % set the active ME to be the Sub ME.
                self.MultiElectrode = self.MultiElectrode.getSubElectrode4ElIdx(channelidx);
                self.activeChannels = self.MultiElectrode.parentElectrodeIndex;
            end
        end
        %------------------------------------------------------------------
        function new = copy(self)
            tmpName = ['temp__' num2str(round(rand*10e10)) '.mat'];
            save(tmpName, 'self');
            Foo = load(tmpName);
            new = Foo.self;
            delete(tmpName);
        end

        %------------------------------------------------------------------
        function setMultiElectrode(self, ME)
            assert(ME.getNElectrodes() == self.getNChannels(), 'The multielectrode must have the same number of channels as me!');
            self.MultiElectrode = ME;
        end
        %------------------------------------------------------------------
        function ME = getMultiElectrode(self)
            ME = self.MultiElectrode;
        end
        %------------------------------------------------------------------
        function b = hasMultiElectrode(self)
            b = ~isempty(self.MultiElectrode);
        end
        %------------------------------------------------------------------
        function wf = getWaveform(self, t, cutLeft, cutLength, channelindex)
            assert(~isempty(t), 't must not be empty!');
            assert( isa(t, 'double'), 'Spiketimes need to be in double format');
            if nargin == 4 || isempty(channelindex)
                channelindex = 1:self.size(2);
            end
            assert(all(channelindex > 0), 'channel index out of bounds!');
            assert(all(channelindex <= self.size(2)), 'channel index out of bounds!');              
            % convert to actual channel indices !
            if ~isempty(self.activeChannels)
                channelindex = self.activeChannels(channelindex);
            end      
          
            if nargin < 6 blockwise = true; end
            
            nT = length(t);
            wf = zeros(nT, length(channelindex)*cutLength);
            
            % Sort the spiketimes and keep the unsorting index:
            [t idx_sort] = sort(t);
            [~, idx_unsort] = sort(idx_sort);
            
            tr = round(t);
            if any(tr~=t)
                warning('t contains non integer values. Subsample shifting while spike cutting is currently not supported. t will be rounded!');
            end
            t = tr;
            
            % Cut away spikes that are at the borders
            while ~isempty(t) && (t(1) <= cutLeft)
                t = t(2:end);
                idx_sort = idx_sort(2:end);
                idx_unsort = idx_unsort(2:end);
            end
            while ~isempty(t) && (t(end)-cutLeft+cutLength+1 > self.size(1))
                t(end) = [];
                idx_sort(end) = [];
                idx_unsort(end) = [];
            end
            
            if isempty(t)
                return
            end
            
            t1 = t-cutLeft;
            t2 = t-cutLeft+cutLength-1;
            nCut = length(t1);
            if nCut < nT
                warning('Could not cut all waveforms, replacing the ones at the edges with zeros (%d from %d cut)', nCut, nT);
            end
            
            wf_t = self.getWaveform_(nCut, channelindex, cutLength, t1, t2);
            % Copy cut waveforms back into original order, leaving those
            % that were too close to the borders zero. Might not the best
            % way to deal with this, but hey, what can we do?            
            wf(idx_sort,:) = wf_t;
        end
        
        %------------------------------------------------------------------
        % Do this in an individual function since certain
        % DataSourceInterface Implementations might want to overwrite the
        % exact way in which the file is accessed. E.g., the H5 matrix,
        % cannot access an irregular index set and has to do this either
        % chunked or individual
        function wf = getWaveform_(self, nCut, channelindex, cutLength, t1, t2)
            wf = zeros(nCut, length(channelindex)*cutLength);
            % Build complete index set and access data with a single
            % entry
            IDX = zeros(1, nCut*cutLength);
            for i = 1:nCut
                idx_idx = (i-1)*cutLength +1 : i*cutLength;
                IDX(idx_idx) = t1(i):t2(i); 
            end 
            X = self.getData_(IDX, channelindex);
            X = reshape(X, [cutLength nCut length(channelindex)]);
            for i = 1:nCut
                wf(i,:) = mysort.wf.m2v(squeeze(X(:, i, :))');
            end
        end
    
        %------------------------------------------------------------------
        function X = getData(self, timeindex, channelindex)
            if nargin < 3
                channelindex = 1:self.size(2);
                if nargin < 2
                    timeindex = 1:self.size(1);
                end
            end
            if ischar(timeindex) && strcmp(timeindex, ':')
                timeindex = 1:self.size(1);
            end
            if ischar(channelindex) && strcmp(channelindex, ':')
                channelindex = 1:self.size(2);
            end
            if ~isempty(self.activeChannels)
                channelindex = self.activeChannels(channelindex);
            end
            X = self.getData_(timeindex, channelindex);
            
            % If this data source should return the residual, remove
            % sorting
            if self.bReturnSortingResiduals && self.hasSpikeSorting()
                S = self.getActiveSpikeSorting();
                if ~isempty(S) && isa(S, 'mysort.spiketrain.SpikeSortingContainer')
                    gdf = S.getGdf(timeindex(1), timeindex(end));
                    if ~isempty(gdf)
                        %                         unitNames = unique(gdf(:,1));
                        T = S.getTemplateWaveforms();
                        %                         T = mysort.wf.pruneTemplates(T, ...
                        %                             'maxChannels', self.maxTemplateChannelsToPlot,...
                        %                             'minChannels', 1, ...
                        %                             'absThreshold', self.templateChannelThreshold,...
                        %                             'setInvalidChannelsTo', nan);
                        cutLeft = S.getTemplateCutLeft();
                        gdf(:,2) = gdf(:,2) - timeindex(1)+1;
                        [x Y] = mysort.wf.templateSpikeTrainData(T, gdf, cutLeft);
                        removeidx = isnan(x) | x < 1 | x > size(X,1);
                        x(removeidx) = [];
                        Y(:,removeidx) = [];
                        X(x,:) = X(x,:) - Y';
                    end
                end
            end
        end
    end
end