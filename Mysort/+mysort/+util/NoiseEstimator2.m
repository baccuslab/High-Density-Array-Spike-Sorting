
classdef NoiseEstimator2 < mysort.util.CacheableClass
    properties
        DH               % DataHandle to manage data access
        noiseEpochs      % matrix containing noise epochs as rows
        xcovs            % cell matrix containing buffered cross 
                         % covariance functions
        noiseEpochBuffer % buffer for the prepared noise epochs that are
                         % shared between to channels. cell matrix
    end
    
    methods
        %%% ----------------CONSTRUCTOR---------------------------     
        function self = NoiseEstimator2(varargin)
            self = self@mysort.util.CacheableClass(varargin{:});
            if self.was_restored; return; end
        end
        
        %%% ----------------CONSTRUCTOR INIT---------------------- 
        function init(self, DH, noiseEpochs, varargin)
            self.P.minCondNumber = 10000;
            self.P.diagonalLoading = 'DL'; % may be 'DL', 'DSL', 'none'
            suppress_warning = 1;
            self.P = mysort.util.parseInputs(self.P, 'NoiseEstimator2', varargin, suppress_warning);
            
            self.DH = DH;
            self.noiseEpochs = noiseEpochs;
            self.xcovs = cell(DH.nC);
            self.noiseEpochBuffer = cell(DH.nC);
            assert(self.P.minCondNumber > 1, 'Condition number cant be lower than 1!');
        end        
        
        %%% ------------------------------------------------------  
        function setNoiseEpochs(self, NE)
            self.noiseEpochs = NE;
        end
        
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
        function Cblock = getBlockNoiseCovarianceMatrix(self, Tf, c1, c2)
            xcov = self.getCrossCovarianceFunction(Tf, c1, c2);
            Cblock = toeplitz( ...
                 xcov(Tf:2*Tf-1), ...
                 xcov(Tf:-1:1) ...
                 );
        end
        
        %%% ------------------------------------------------------
        function xcov = getCrossCovarianceFunction(self, Tf, c1, c2)
            transpose = 0;
            if c1 > c2
                tmp = c1;
                c1 = c2;
                c2 = tmp;
                transpose = 1;
            end
            
            % Load buffer if possible
            if isempty(self.xcovs{c1, c2}) ||...
               ~isfield(self.xcovs{c1,c2}, 'Tf') || ...
               self.xcovs{c1,c2}.Tf < Tf
                xcov = self.calcXCov(Tf, c1, c2);
                self.xcovs{c1, c2}.xcov = xcov;
                self.xcovs{c1, c2}.Tf = Tf;
            else
                middleIdx = (length(self.xcovs{c1, c2}.xcov)-1)/2 +1;
                xcov = self.xcovs{c1, c2}.xcov( middleIdx - (Tf-1):middleIdx + (Tf-1) );
            end
            
            if transpose == 1
                xcov = xcov(end:-1:1);
            end
        end

        %%% ------------------------------------------------------
        function xcov = calcXCov(self, Tf, c1, c2)
            [noiseEpochs epochLength] = self.prepareNoiseEpochs(Tf, c1, c2);
    
            totalNoiseEpochLength = sum(epochLength);
            xcov = zeros(1, 2*(Tf-1) +1 );
            if c1 == c2
                for i=1:size(noiseEpochs,1)
                    xcov = xcov +  xcorr(...
                        self.DH.getData('start', noiseEpochs(i,1), ...
                                        'stopp', noiseEpochs(i,2), ...
                                        'channels', c1), ...
                        Tf-1, 'none');
                end
            else
                for i=1:size(noiseEpochs,1)
                    X = self.DH.getData('start', noiseEpochs(i,1), ...
                                        'stopp', noiseEpochs(i,2), ...
                                        'channels', [c1 c2]);
                    xcov = xcov + xcorr(...
                        X(1,:), X(2,:), Tf-1, 'none');
                end
            end
            xcov = xcov/totalNoiseEpochLength;            
        end
        
        %%% ------------------------------------------------------
        function [noiseEpochs epochLength] = prepareNoiseEpochs(self, Tf, c1, c2)
            % TODO: Allow for channel wise noise epochs and use c1 and c2
            % to get the common noise epochs of those two channels
            if ~isempty(self.noiseEpochBuffer{c1, c2}) || ...
               ~isfield(self.noiseEpochBuffer{c1, c2}, 'Tf') || ...
               self.noiseEpochBuffer{c1, c2}.Tf < Tf
                epochLength = self.noiseEpochs(:,2)-self.noiseEpochs(:,1)+1;
                noiseEpochs = self.noiseEpochs(epochLength > 3*Tf,:);
                epochLength = epochLength( epochLength > 3*Tf);
                self.noiseEpochBuffer{c1, c2}.Tf = Tf;
                self.noiseEpochBuffer{c1, c2}.noiseEpochs = noiseEpochs;
                self.noiseEpochBuffer{c1, c2}.epochLength = epochLength;
            else
                noiseEpochs = self.noiseEpochBuffer{c1, c2}.noiseEpochs;
                epochLength = self.noiseEpochBuffer{c1, c2}.epochLength;
            end            
        end
        
        %%% ------------------------------------------------------
        function c = cond(self, Tf, channels)
           d = self.eig(Tf, channels);
           if any(d == 0)   % Handle singular matrix
               c = Inf;
           else
               c = max(d)./min(d);
           end
        end
        
        %%% ------------------------------------------------------
        function [d_V, D] = eig(self, Tf, channels)
            if nargout == 1
                d_V = eig(self.getNoiseCovarianceMatrix(Tf, channels));
            else
                [d_V, D] = eig(self.getNoiseCovarianceMatrix(Tf, channels));
            end
        end
    end
end