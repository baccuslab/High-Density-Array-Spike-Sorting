classdef Covest < mysort.util.DebuggableClass
    properties
        dataSource
        xcovs
        CCol
    end
    
    methods
        %------------------------------------------------------------------
        function self = Covest(X, varargin)
            self = self@mysort.util.DebuggableClass(varargin{:});
            self.P.maxDist = 50; % micro meter
            self.P.maxLag  = 20; % samples, should be = Tf-1
            self.P.noiseEpochs = [];
            self.P.delayComputation = 0;
            self.P.forceMethod = []; % possible: "xcorr", "matmul"
            self.P = mysort.util.parseInputs(self.P, '', varargin);
            self.dataSource = mysort.datasource.check(X);
            
            if isempty(self.P.noiseEpochs)
                self.P.noiseEpochs = [1 self.dataSource.getLen()];
            end
            
            if ~self.P.delayComputation
                self.xcovs = self.calcXCovs();
                self.CCol  = self.buildCColumn();
            end
        end
        
        %------------------------------------------------------------------
        function xcovs = calcXCovs(self)
            if strcmp(self.P.forceMethod, 'matmul')
                xcovs = self.calcXCovsWithMatMul();
            elseif strcmp(self.P.forceMethod, 'xcorr');
                xcovs = self.calcXCovsWithXCorr();
            elseif isempty(self.P.forceMethod)
                if self.P.maxLag <= 30 % this is the break-even point. see 
                                       % \mysortpackage\tests\tests\noise\covestTest2.m
                    xcovs = self.calcXCovsWithMatMul();
                else 
                    xcovs = self.calcXCovsWithXCorr();
                end                
            else
                error('Method not defined: %s', self.P.forceMethod);
            end
        end
        %% FOR LONG time lags
        %------------------------------------------------------------------
        function xcovs = calcXCovsWithXCorr(self)
            nC = self.dataSource.getNChannel();
            xcovs = cell(nC, nC);
            for i=1:nC
                for j=i:nC
                    if isa(self.dataSource, 'mysort.datasource.MultiElectrodeInterface')
                        if self.dataSource.getDistance(i,j) < self.P.maxDist
                            xcovs{i,j} = self.calcXCovBetweenChannelWithXCorr(i,j);
                        end                            
                    else
                        xcovs{i,j} = self.calcXCovBetweenChannelWithXCorr(i,j);
                    end
                end
            end        
        end
        
        %------------------------------------------------------------------
        function xcov = calcXCovBetweenChannelWithXCorr(self, c1, c2)
            maxLag = self.P.maxLag;
            xcov = zeros(1, 2*maxLag +1 );
            totalNoiseEpochLength = sum(mysort.epoch.length(self.P.noiseEpochs));
            if c1 == c2
                for i=1:size(self.P.noiseEpochs,1)
                    xcov = xcov +  xcorr(...
                        self.dataSource.getData(...
                            [self.P.noiseEpochs(i,1) self.P.noiseEpochs(i,2)], ...
                            c1), ...
                        maxLag, 'none');
                end
            else
                for i=1:size(self.P.noiseEpochs,1)
                    X = self.dataSource.getData(...
                        [self.P.noiseEpochs(i,1) self.P.noiseEpochs(i,2)], ...
                        [c1 c2]);
                    xcov = xcov + xcorr(...
                        X(1,:), X(2,:), maxLag, 'none');
                end
            end
            xcov = xcov/totalNoiseEpochLength;              
        end
        
        %% FOR SMALL timelags
        %------------------------------------------------------------------
        function xcovs = calcXCovsWithMatMul(self)
            nC = self.dataSource.getNChannel();
            xcovs = cell(nC, nC);
            maxlag = self.P.maxLag;
            totalNoiseEpochLength = sum(mysort.epoch.length(self.P.noiseEpochs));
            % do this only once since "isa" and getDistance are very slow 
            % when called often (50% of total computation time)
            bIsMea = isa(self.dataSource, 'mysort.datasource.MultiElectrodeInterface');
            bufferedChannelPairs = prepareChannelPairs();
            X = self.dataSource.getData(self.P.noiseEpochs);
            EL = mysort.epoch.length(self.P.noiseEpochs);
            ELcumsum = [0; cumsum(EL)];
            for ne=1:size(self.P.noiseEpochs,1)
                for c1=1:nC
                    xcovs{c1,c1} = calcXCovBetweenChannelMatMul(...
                        X(c1, ELcumsum(ne)+1:ELcumsum(ne+1)), ...
                        X(c1, ELcumsum(ne)+1:ELcumsum(ne+1)));
                    for c2=c1+1:nC
                        if bufferedChannelPairs(c1,c2)
                            if isempty(xcovs{c1,c2})
                                xcovs{c1,c2} = calcXCovBetweenChannelMatMul(...
                                    X(c1, ELcumsum(ne)+1:ELcumsum(ne+1)),...
                                    X(c2, ELcumsum(ne)+1:ELcumsum(ne+1)));
                            else
                                xcovs{c1,c2} = xcovs{c1,c2} + ...
                                    calcXCovBetweenChannelMatMul(...
                                    X(c1, ELcumsum(ne)+1:ELcumsum(ne+1)),...
                                    X(c2, ELcumsum(ne)+1:ELcumsum(ne+1)));
                            end
                        end
                    end
                end
            end
            % Normalize covariances
            for c1 = 1:nC
                for c2 = c1:nC
                    xcovs{c1,c2} = xcovs{c1,c2}/totalNoiseEpochLength;
                end
            end
            %----
            function xcov = calcXCovBetweenChannelMatMul(d1, d2)
                xcov = zeros(1, 2*maxlag+1);
                tau0 = maxlag+1;
                xcov(tau0) = d1*d2';
                for tau = 1:maxlag
                    xcov(tau0 + tau) = d1(tau+1:end)*d2(1:end-tau)';
                    xcov(tau0 - tau) = d1(1:end-tau)*d2(tau+1:end)';
                end
            end
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
        
        %% Rest        
        %------------------------------------------------------------------
        function ccol = buildCColumn(self, maxLag)
            if ~exist('maxLag', 'var')
                maxLag = self.P.maxLag;
            end
            ccol = mysort.noise.xcov2ccol(self.xcovs, maxLag);
        end       
        
        %------------------------------------------------------------------
        function x = invMul(self, y)
            % solves C*x = y; <=> x = inv(C)*y
            % if y is a matrix invMul operates on the ROWS of y
            if size(y,1) > 1 && size(y,2) > 1
                % y is a matrix, solve every row individually
                x = zeros(size(y));
                for i=1:size(y,1)
                    x(i,:) = self.invMul(y(i,:));
                end
            else
                % y is a vector, keep orientation
                turn = 0;
                if size(y,1) == 1
                    turn = 1;
                    y = y';
                end
                nC = self.dataSource.getNChannel();
                maxLag = length(y)/nC -1;
                ccol = self.buildCColumn(maxLag);
                x = matlabfilecentral.block_levinson(y, ccol);

                if turn
                    x = x';
                end
            end
        end
        
        %% Get ChannelEmbedding
        %%% ------------------------------------------------------  
        function C = getNoiseCovarianceMatrix(self, Tf, channels)
            C = zeros(Tf*length(channels), Tf*length(channels));
            for row=1:length(channels)
                row_start_idx = (row-1)*Tf +1;
                row_stopp_idx  = row*Tf;
                for column=1:length(channels)
                    col_start_idx = (column-1)*Tf +1;
                    col_stopp_idx  = column*Tf;
                    
                    C(row_start_idx:row_stopp_idx, col_start_idx:col_stopp_idx) = ...
                        self.getBlockNoiseCovarianceMatrix(Tf, channels(row), channels(column));
                end
            end
        end
        
        %%% ------------------------------------------------------
        function C = getNoiseCovarianceMatrixTimeEmbed(self, Tf, channels)
            C = zeros(Tf*length(channels), Tf*length(channels));
            TfInXcov = (length(self.xcovs{1,1})-1)/2;
            d = TfInXcov-Tf;

            for row=1:length(channels)
                row_start_idx = (row-1)*Tf +1;
                row_stopp_idx  = row*Tf;
                for column=row:length(channels)
                    col_start_idx = (column-1)*Tf +1;
                    col_stopp_idx  = column*Tf;
                    xcov = self.xcovs{row, column};
                    if isempty(xcov)
                        xcov = zeros(1, 2*TfInXcov+1);
                    end

                    if d>=0
                        xcov = xcov(d+1:end-d);
                    else
                        xcov = [zeros(1,-d) xcov zeros(1,-d)];
                    end
                    mIdx = (length(xcov)+1)/2;
                    Cblock = toeplitz( ...
                             xcov(mIdx:mIdx+Tf-1), ...
                             xcov(mIdx:-1:mIdx-Tf+1) ...
                             );
                    C(row_start_idx:row_stopp_idx, col_start_idx:col_stopp_idx) = ...
                        Cblock;
                    if column>row
                        C(col_start_idx:col_stopp_idx, row_start_idx:row_stopp_idx) = ...
                                    Cblock';
                    end
                end
            end
        end        
    end
end