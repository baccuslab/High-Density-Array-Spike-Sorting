classdef FakeCovest < mysort.noise.Covest
    properties
        iC
        bInvertWholeMatrix
    end
    
    methods
        %------------------------------------------------------------------
        function self = FakeCovest(X, bInvertWholeMatrix, varargin)
            self = self@mysort.noise.Covest(X, 'delayComputation', 1, varargin{:});
            self.bInvertWholeMatrix = bInvertWholeMatrix;
            self.xcovs = self.calcXCovs();
        end
        
        %------------------------------------------------------------------
        function xcovs = calcXCovs(self)
            nC = self.dataSource.getNChannel();
            xcovs = zeros(nC, nC);
            maxlag = 0;
            totalNoiseEpochLength = sum(mysort.epoch.length(self.P.noiseEpochs));
            % do this only once since "isa" and getDistance are very slow 
            % when called often (50% of total computation time)
            bIsMea = isa(self.dataSource, 'mysort.datasource.MultiElectrodeInterface');
%             bufferedChannelPairs = prepareChannelPairs();
            X = self.dataSource.getData(self.P.noiseEpochs);
            disp('Computing Cov...'); tic
            xcovs = cov(X');
            toc
            disp('Inverting Cov...'); tic
            if self.bInvertWholeMatrix
                self.iC = inv(xcovs);
            else
                self.iC = diag(diag(xcovs).^-1);
            end
            toc
            
%             for c1=1:nC
%                 xcovs(c1,c1) = X(c1,:)*X(c1,:)' / totalNoiseEpochLength;
%                 for c2=c1+1:nC
%                     if bufferedChannelPairs(c1,c2)
%                         xcovs(c1,c2) = X(c1,:)*X(c2,:)' / totalNoiseEpochLength;
%                     end
%                 end
%             end
            
            %----
            function CP = prepareChannelPairs()
                CP = ones(nC,nC);
                if ~bIsMea
                    return
                end
                % channel 2 data is only needed, if the
                % distance between c1 and c2 is close enough.
                % if we dont know that distance, we also have 
                % to calculate the xcov
                for cc1=1:nC
                    for cc2=cc1:nC
                        CP(cc1,cc2) = self.dataSource.getDistance(cc1,cc2) < self.P.maxDist;
                    end
                end
            end             
        end
        %% FOR LONG time lags
        %------------------------------------------------------------------
        function xcovs = calcXCovsWithXCorr(self)
            error('Dont use this function');
        end
        
        %------------------------------------------------------------------
        function xcov = calcXCovBetweenChannelWithXCorr(self, c1, c2)
            error('Dont use this function');            
        end
        
        %% FOR SMALL timelags
        %------------------------------------------------------------------
        function xcovs = calcXCovsWithMatMul(self)
            error('Dont use this function');
        end
        
        %% Rest        
        %------------------------------------------------------------------
        function ccol = buildCColumn(self, maxLag)
            if ~exist('maxLag', 'var')
                maxLag = self.P.maxLag;
            end
            
        end       
        
        %------------------------------------------------------------------
        function x = invMul(self, y)
        	nC = self.dataSource.getNChannel();
            nT = size(y,1);
            x = zeros(size(y));
            for i=1:nT
                x(i,:) = mysort.wf.m2v(self.iC*mysort.wf.v2m(y(i,:), nC), nC);
            end
        end
        
        %% Get ChannelEmbedding
        %%% ------------------------------------------------------  
        function C = getNoiseCovarianceMatrix(self, Tf, channels)
            error('Dont use this function');
        end
        
        %%% ------------------------------------------------------
        function C = getNoiseCovarianceMatrixTimeEmbed(self, Tf, channels)
            error('Dont use this function');
        end        
    end
end