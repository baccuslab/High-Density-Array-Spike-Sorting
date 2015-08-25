classdef Covest2 < handle
    properties
        DS
        CCol
        P
    end
    
    methods
        %------------------------------------------------------------------
        function self = Covest2(X, varargin)
            self.P.maxDist = []; % micro meter
            self.P.maxLag  = 79; % samples, should be = Tf-1
            self.P.maxSamples = 100000;
            self.P.maxSamplesPerEpoch = 100000;
            self.P.isCCol = 0;
            self.P.noiseEpochs = [];
            self.P.forceMethod = 'xcorr'; % possible: 'xcorr', 'matmul'
            self.P = mysort.util.parseInputs(self.P, varargin, 'error');
            if isstruct(X)
                self.fromStruct(X);
                return
            end
            
            % remove short
            self.P.noiseEpochs = mysort.epoch.removeShort(self.P.noiseEpochs, self.P.maxLag);
            
            if isempty(self.P.noiseEpochs)
                self.P.noiseEpochs = [1 min(self.P.maxSamples, size(X,1))];
            end
            
            NoiseEpochLens = mysort.epoch.length(self.P.noiseEpochs);
            csum = cumsum(NoiseEpochLens);
            lastidx = find(csum>self.P.maxSamples, 1);
            if isempty(lastidx)
               lastidx = size(self.P.noiseEpochs,1);
            end
            self.P.noiseEpochs = self.P.noiseEpochs(1:lastidx,:);
            
            % Make sure there is not one very long epoch
            self.P.noiseEpochs = mysort.epoch.makeMaxLen(self.P.noiseEpochs, self.P.maxSamplesPerEpoch);
            
            if self.P.isCCol
                self.CCol = X;
                nC = size(self.CCol,2);
                self.P.maxLag = min(size(X,1)/nC, self.P.maxLag);
            elseif ~isempty(X)
%                 self.DS = mysort.ds.Matrix(X(1:self.P.noiseEpochs(end,2),:), X.getSampleRate(), X.MultiElectrode);
                self.DS = X;
                xcovs = self.calcXCovs();
                self.CCol = mysort.noise.xcov2ccol(xcovs, self.P.maxLag);
            end
        end
        %------------------------------------------------------------------
        function S = toStruct(self)
            S.P = self.P;
            S.CCol = self.CCol;
        end
        %------------------------------------------------------------------
        function fromStruct(self, S)
            assert(~isempty(S.CCol), 'CCol must not be empty!');
            self.P = S.P;
            self.CCol = S.CCol;
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
            nC = size(self.DS,2);
            xcovs = cell(nC, nC);
            maxLag = self.P.maxLag;
            totalNoiseEpochLength = sum(mysort.epoch.length(self.P.noiseEpochs));
            CP = ones(nC,nC);
            for c1 = 1:nC
                for c2 = c1:nC
                    xcovs{c1, c2} = zeros(1, 2*maxLag +1 );
                    if ~isempty(self.P.maxDist)
                        CP(c1,c2) = self.DS.MultiElectrode.getDistance(c1,c2) < self.P.maxDist;
                    end
                end
            end
            % only load data of one chunk once!
            for e=1:size(self.P.noiseEpochs,1)
                X = self.DS(self.P.noiseEpochs(e,1):self.P.noiseEpochs(e,2), :);
                for c1=1:nC
                    for c2=c1:nC
                        if c1==c2
                            xcovs{c1,c1} = xcovs{c1,c1} + xcorr(X(:, c1), maxLag, 'none')';
                        elseif CP(c1,c2)
                            xcovs{c1,c2} = xcovs{c1,c2} + xcorr(X(:,c1), X(:,c2), maxLag, 'none')';
                        end
                    end
                end        
            end
            for c1 = 1:nC
                for c2 = c1:nC
                    xcovs{c1, c2} = xcovs{c1, c2}/totalNoiseEpochLength; 
                end
            end            
        end
        
        %% FOR SMALL timelags
        %------------------------------------------------------------------
        function xcovs = calcXCovsWithMatMul(self)
            error('Dont use this code, it contains a bug!');
            nC = size(self.DS,2);
            xcovs = cell(nC, nC);
            maxlag = self.P.maxLag;
            NEL = mysort.epoch.length(self.P.noiseEpochs);
            totalNoiseEpochLength = sum(NEL);
            % do this only once since "isa" and getDistance are very slow 
            % when called often (50% of total computation time)
            bufferedChannelPairs = prepareChannelPairs();
            X = zeros(totalNoiseEpochLength, nC);
            s1 = 1;
            for i=1:size(self.P.noiseEpochs,1)
                X(s1:s1+NEL(i)-1,:) = self.DS(self.P.noiseEpochs(i,1):self.P.noiseEpochs(i,2),:);
                s1 = s1+NEL(i);
            end
            EL = mysort.epoch.length(self.P.noiseEpochs);
            ELcumsum = [0; cumsum(EL)];
            for ne=1:size(self.P.noiseEpochs,1)
                for c1=1:nC
                    xcovs{c1,c1} = calcXCovBetweenChannelMatMul(...
                        X(ELcumsum(ne)+1:ELcumsum(ne+1), c1), ...
                        X(ELcumsum(ne)+1:ELcumsum(ne+1), c1));
                    for c2=c1+1:nC
                        if bufferedChannelPairs(c1,c2)
                            if isempty(xcovs{c1,c2})
                                xcovs{c1,c2} = calcXCovBetweenChannelMatMul(...
                                    X(ELcumsum(ne)+1:ELcumsum(ne+1), c1),...
                                    X(ELcumsum(ne)+1:ELcumsum(ne+1), c2));
                            else
                                xcovs{c1,c2} = xcovs{c1,c2} + ...
                                    calcXCovBetweenChannelMatMul(...
                                    X(ELcumsum(ne)+1:ELcumsum(ne+1), c1),...
                                    X(ELcumsum(ne)+1:ELcumsum(ne+1), c2));
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
                xcov(tau0) = d1'*d2;
                for tau = 1:maxlag
                    xcov(tau0 + tau) = d1(tau+1:end)'*d2(1:end-tau);
                    xcov(tau0 - tau) = d1(1:end-tau)'*d2(tau+1:end);
                end
            end
            %----
            function CP = prepareChannelPairs()
                CP = ones(nC,nC);
                % channel 2 data is only needed, if the
                % distance between c1 and c2 is close enough.
                % if we dont know that distance, we also have 
                % to calculate the xcov
                for cc1=1:nC
                    CP(cc1,cc1) = true;
                    for cc2=cc1+1:nC
                        CP(cc1,cc2) = self.DS.MultiElectrode.getDistance(cc1,cc2) < self.P.maxDist;
                    end
                end
            end
        end
        
        %% Rest        
                
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
                nC = self.DS.getNChannels();
                maxLag = length(y)/nC -1;
                ccol = self.CCol(1:nC*(maxLag+1),:);
                x = matlabfilecentral.block_levinson(y, ccol);

                if turn
                    x = x';
                end
            end
        end
    end
end