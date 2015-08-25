
classdef NoiseEstimator < mysort.util.DebuggableClass
    properties
        DH               % DataHandle to manage data access
        noiseEpochs      % matrix containing noise epochs as rows
        xcovs            % cell matrix containing buffered cross 
                         % covariance functions
        noiseEpochBuffer % buffer for the prepared noise epochs that are
                         % shared between to channels. cell matrix
        eigenDecompositionHashD  % Java Hashtable to buffer the eigen
                                 % decompositions of subcovariance matrices
        eigenDecompositionHashV
        covarianceMatrixHash
        full_noise_cov_matrix    % alternative way to use this class is to 
                                 % provide the full noise covariance matrix
                                 % if available
        Tf_of_full_cov_matrix    % only used if full noise cov available
        nC               % number of channels
    end
    
    methods
        %%% ----------------CONSTRUCTOR---------------------------     
        function self = NoiseEstimator(DH_or_C_or_xcovs, noiseEpochs_or_Tf, varargin)
            self = self@mysort.util.DebuggableClass(varargin{:});
        
            % TODO: KANÄLE SORTIEREN VOR DEM PUFFER !!!! 
        
            self.P.minCondNumber = 10000;
            self.P.diagonalLoading = 'DL'; % may be 'DL', 'DSL', 'none'
            suppress_warning = 1;
            self.P = mysort.util.parseInputs(self.P, 'NoiseEstimator', varargin, suppress_warning);
            self.DH = [];
            self.noiseEpochs = [];      
            self.noiseEpochBuffer = [];
            
            if iscell(DH_or_C_or_xcovs)
                % DH_or_C_or_xcovs is a xcov cell
                % noiseEpochs_or_Tf is not defined
                self.xcovs = DH_or_C_or_xcovs;
                self.nC = size(DH_or_C_or_xcovs,1);
            elseif isnumeric(DH_or_C_or_xcovs) && isnumeric(noiseEpochs_or_Tf)
                % DH_or_C_or_xcovs is a Covariance Matrix
                % noiseEpochs_or_Tf is Tf
                self.nC = size(DH_or_C_or_xcovs,1)/noiseEpochs_or_Tf;
                assert(round(self.nC) == self.nC, 'Tf must divide size of full covariance matric!');
                assert(size(DH_or_C_or_xcovs,1) == size(DH_or_C_or_xcovs,2), 'Covariance matrix must be square!');
                self.full_noise_cov_matrix = DH_or_C_or_xcovs;
                self.Tf_of_full_cov_matrix = noiseEpochs_or_Tf;
                self.xcovs = [];
            elseif isa(DH_or_C_or_xcovs, 'mysort.datafile.DataFileInterface')
                % DH_or_C_or_xcovs is a DataHandle
                % noiseEpochs_or_Tf are noise epochs               
                self.DH = DH_or_C_or_xcovs;
                self.nC = self.DH.nC;
                self.noiseEpochs = noiseEpochs_or_Tf;
                self.xcovs = cell(self.DH.nC);
                self.noiseEpochBuffer = cell(self.DH.nC);
            else
                error('not implemented')
            end
            assert(self.P.minCondNumber > 1, 'Condition number cant be lower than 1!');
            self.eigenDecompositionHashD = java.util.Hashtable;
            self.eigenDecompositionHashV = java.util.Hashtable;  
            self.covarianceMatrixHash    = mysort.util.HashBufferedFunction(@(Tf, chan) self.getNoiseCovarianceMatrix_(Tf, chan));
        end        
        
        %%% ------------------------------------------------------  
        function setNoiseEpochs(self, NE)
            self.noiseEpochs = NE;
        end
        
        %%% ------------------------------------------------------  
        function C = getNoiseCovarianceMatrix(self, Tf, channels)
            if ~exist('channels', 'var') || isempty(channels); channels = 1:self.nC; end
%             hash = [num2str(Tf) '_' num2str(channels)];
%             hash = int64(java.lang.String(hash).hashCode);
%             C = self.covarianceMatrixHash.call(hash, Tf, channels);
            C = self.getNoiseCovarianceMatrix_(Tf, channels);
        end
        
        %%% ------------------------------------------------------  
        function C = getNoiseCovarianceMatrix_(self, Tf, channels)
            if ~exist('channels', 'var') || isempty(channels); channels = 1:self.nC; end
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
            if isempty(self.full_noise_cov_matrix)
                xcov = self.getCrossCovarianceFunction(Tf, c1, c2);
                Cblock = toeplitz( ...
                     xcov(Tf:2*Tf-1), ...
                     xcov(Tf:-1:1) ...
                     );
            else
                % the full noise covariance matrix is available. use it
                assert(Tf <= self.Tf_of_full_cov_matrix, sprintf('The full noise covariance matrix was provided with a Tf of %d. You cant compute a bigger one (Tf=%d)\n', self.Tf_of_full_cov_matrix, Tf));
                s1 = (c1-1)*self.Tf_of_full_cov_matrix+1;
                s2 = (c2-1)*self.Tf_of_full_cov_matrix+1;
                Cblock = self.full_noise_cov_matrix(s1:s1+Tf-1, s2:s2+Tf-1);
            end
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
            assert(~isempty(self.noiseEpochBuffer), 'noiseEpochBuffer was not correctly initialized!');
            if ~isempty(self.noiseEpochBuffer{c1, c2}) || ...
               ~isfield(self.noiseEpochBuffer{c1, c2}, 'Tf') || ...
               self.noiseEpochBuffer{c1, c2}.Tf < Tf
                epochLength = self.noiseEpochs(:,2)-self.noiseEpochs(:,1)+1;
                noiseEpochs = self.noiseEpochs(epochLength >= Tf,:);
                epochLength = epochLength( epochLength >= Tf);
                self.noiseEpochBuffer{c1, c2}.Tf = Tf;
                self.noiseEpochBuffer{c1, c2}.noiseEpochs = noiseEpochs;
                self.noiseEpochBuffer{c1, c2}.epochLength = epochLength;
            else
                noiseEpochs = self.noiseEpochBuffer{c1, c2}.noiseEpochs;
                epochLength = self.noiseEpochBuffer{c1, c2}.epochLength;
            end            
        end
        
        %%% ------------------------------------------------------
        function c = cond(self, varargin)
           d = self.eig(varargin{:});
           if any(d == 0)   % Handle singular matrix
               c = Inf;
           else
               c = max(d)./min(d);
           end
        end
        
        %%% ------------------------------------------------------
        function [d_V, D] = eig(self, Tf, channels)
            if ~exist('channels', 'var') || isempty(channels); channels = 1:self.nC; end
            if nargout < 2
                if self.eigenDecompositionHashD.containsKey(num2str([Tf channels]))
                    d_V = self.eigenDecompositionHashD.get(num2str([Tf channels]));
                else
                    d_V = eig(self.getNoiseCovarianceMatrix(Tf, channels));
                    assert(~any(imag(d_V)), 'Covariance matrix has imaginary eigenvalues!');
                    assert(~any(d_V<0), 'Covariance matrix has negative eigenvalues!');
                    self.eigenDecompositionHashD.put(num2str([Tf channels]), d_V);
                end                
            else
                if self.eigenDecompositionHashV.containsKey(num2str([Tf channels]))
                    D = diag(self.eigenDecompositionHashD.get(num2str([Tf channels])));
                    d_V = self.eigenDecompositionHashV.get(num2str([Tf channels]));
                else
                    [d_V, D] = eig(self.getNoiseCovarianceMatrix(Tf, channels));
                    assert(~any(imag(diag(D))), 'Covariance matrix has imaginary eigenvalues!');                    
                    assert(~any(diag(D)<0), 'Covariance matrix has negative eigenvalues!');
                    self.eigenDecompositionHashD.put(num2str([Tf channels]), diag(D));
                    self.eigenDecompositionHashV.put(num2str([Tf channels]), d_V);
                end                
            end
        end
        %%% ------------------------------------------------------
        function Cinv = inv(self, varargin)
            [V, D] = self.eig(varargin{:});
            Cinv = V*diag((1./diag(D)))*V';
        end
        %%% ------------------------------------------------------
        function ch = sqrt(self, varargin)
            [V, D] = self.eig(varargin{:});
            ch = V*(D^(1/2))*V';
        end
        %%% ------------------------------------------------------
        function U = chol(self, varargin)
            U = chol(self.getNoiseCovarianceMatrix(varargin{:}));
        end
        %%% ------------------------------------------------------
        function iU = getPrewhiteningOperator(self, varargin)
            iU = inv(self.chol(varargin{:}));
        end
    end
end