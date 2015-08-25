classdef FilteredDataSourceInterface < mysort.ds.DataSourceInterface
    properties
        filterFactory
        filterObject
        filterObjectInitializationParameters
        useFilter
    end
    
    methods(Abstract)
        getUnfilteredData(self, idx1, idx2)
    end
    
    methods
        %------------------------------------------------------------------
        function self = FilteredDataSourceInterface(filterFactory, useFilter, varargin)
            self = self@mysort.ds.DataSourceInterface(varargin{:});
            self.filterFactory = filterFactory;
            self.useFilter = useFilter;
        end
        %------------------------------------------------------------------
        function X = getData_(self, timeIndex, channelIndex)
            if self.useFilter
                X = self.getFilteredData(timeIndex, channelIndex);
            else
                X = self.getUnfilteredData(timeIndex, channelIndex);
            end
        end
        %------------------------------------------------------------------
        function X = getFilteredData(self, timeIndex, channelIndex, varargin)
            P.filterChunkSize = 500000;
            P.chunkOverlap = 100; % burn-in time for filter
            P.filterProgressDisplay = 'none'; % or 'console' or 'progressbar';
            P.loadProgressDisplay = 'none'; % or 'console' or 'progressbar';
            P = mysort.util.parseInputs(P, varargin, 'error');
            
            t1 = max(1, min(timeIndex)-P.chunkOverlap);            
            offset = min(timeIndex)-t1;
            t2 = min(size(self,1), max(timeIndex)+P.chunkOverlap);
            endOverlap = t2-max(timeIndex);
            X = self.getUnfilteredData(t1:t2, channelIndex, ...
                'progressDisplay', P.loadProgressDisplay);

            
            if 0 %length(timeIndex) > 1000 && matlabpool('size') > 1 && size(X,2)>1
                cps = {};
                for c=1:size(X,2)
                    cps{c} = self.filterFactory.getFilterCopy();
                end
                % use parallel processing if available
                % no, this does not save time!
                parfor c = 1:size(X,2)
                    X(:,c) = filtfilthd(cps{c}, X(:,c));
                end
            else
                cp = self.filterFactory.getFilterCopy();
                X = filtfilthd(cp, X);
            end
            X = X(offset+1:end-endOverlap, :);
            X = X(timeIndex-min(timeIndex)+1,:);
        end        
        %------------------------------------------------------------------
        function initializeFilter(self, electrodeIndex)
            if nargin==1 || isempty(electrodeIndex)
                electrodeIndex = 1:self.size(2);
            end
            if isempty(self.filterObject)
                % The filter is needed. Design now.
                self.filterObject = self.filterFactory.getFilterCopy();
            else
                reset(self.filterObject);
            end
            if isempty(self.filterObjectInitializationParameters)
                fprintf('Initializing filter...');
                nC = length(electrodeIndex);
                startSamples = self.getUnfilteredData(1:min(self.size(1), 10000), 1:nC);
                self.filterObjectInitializationParameters = mean(startSamples, 1);   %first 50ms            
                Y = filter(self.filterObject, startSamples);
                self.filterObjectInitializationParameters = self.filterObject.States;
                fprintf(' done.\n');
            end
            self.filterObject.PersistentMemory=1;
%             m = self.filterObjectInitializationParameters(electrodeIndex);
%             self.filterObject.States = m([1 1],:); %2nd order bandpass (FIXME for others...)
            self.filterObject.States = self.filterObjectInitializationParameters;
        end
    end
end