classdef Covestrev21918 < mysort.util.DebuggableClass
    properties
        dataSource
        xcovs
        CCol
    end
    
    methods
        %------------------------------------------------------------------
        function self = Covestrev21918(X, varargin)
            self = self@mysort.util.DebuggableClass(varargin{:});
            self.P.maxDist = 50; % micro meter
            self.P.maxLag  = 50; % samples
            self.P.noiseEpochs = [];
            self.P = mysort.util.parseInputs(self.P, '', varargin);
            if isstruct(X)
                self.fromStruct(X);
                return
            end
            if isa(X, 'mysort.datasource.DataSourceInterface')
                self.dataSource = X;
            else
                self.dataSource = mysort.datasource.DataSource(X);
            end
            
            if isempty(self.P.noiseEpochs)
                self.P.noiseEpochs = [1 self.dataSource.getLen()];
            end
            
            self.xcovs = self.calcXCovs();
            self.CCol  = self.buildCColumn();
        end
        %------------------------------------------------------------------
        function S = toStruct(self)
            S.P = self.P;
            S.xcovs = self.xcovs;
            S.CCol = self.CCol;
        end
        %------------------------------------------------------------------
        function fromStruct(self, S)
            self.P = S.P;
            self.xcovs = S.xcovs;
            self.CCol = S.CCol;
        end
        
        %------------------------------------------------------------------
        function xcovs = calcXCovs(self)
            nC = self.dataSource.getNChannel();
            xcovs = cell(nC, nC);
            for i=1:nC
                for j=i:nC
                    if isa(self.dataSource, 'mysort.datasource.MultiElectrodeInterface')
                        if self.dataSource.getDistance(i,j) < self.P.maxDist
                            xcovs{i,j} = self.calcXCovBetweenChannel(i,j);
                        end                            
                    else
                        xcovs{i,j} = self.calcXCovBetweenChannel(i,j);
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        function xcov = calcXCovBetweenChannel(self, c1, c2)
            Tf = self.P.maxLag;
            xcov = zeros(1, 2*(Tf-1) +1 );
            totalNoiseEpochLength = sum(self.P.noiseEpochs(:,2)-self.P.noiseEpochs(:,1));
            if c1 == c2
                for i=1:size(self.P.noiseEpochs,1)
                    xcov = xcov +  xcorr(...
                        self.dataSource.getData(self.P.noiseEpochs(i,:), c1), ...
                        Tf-1, 'none');
                end
            else
                for i=1:size(self.P.noiseEpochs,1)
                    X = self.dataSource.getData(self.P.noiseEpochs(i,:), [c1 c2]);
                    xcov = xcov + xcorr(...
                        X(1,:), X(2,:), Tf-1, 'none');
                end
            end
            xcov = xcov/totalNoiseEpochLength;              
        end
        %------------------------------------------------------------------
        function ccol = buildCColumn(self)
            ccol = mysort.noise.xcov2ccol(self.xcovs);
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
                x = matlabfilecentral.block_levinson(y, self.CCol);
                if turn
                    x = x';
                end
            end
        end
    end
end