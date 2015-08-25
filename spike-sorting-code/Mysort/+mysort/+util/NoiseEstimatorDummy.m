
classdef NoiseEstimatorDummy < mysort.util.NoiseEstimator
    properties

    end
    
    methods
        %%% ----------------CONSTRUCTOR---------------------------     
        function self = NoiseEstimatorDummy(varargin)
            self = self@mysort.util.NoiseEstimator(varargin{:});
        end        
    
%         %%% ------------------------------------------------------
%         function [d_V, D] = eig(self, Tf, channels)
%             if ~exist('channels', 'var') || isempty(channels); channels = 1:self.nC; end
%             if nargout < 2
%                 if self.eigenDecompositionHashD.containsKey(num2str([Tf channels]))
%                     d_V = self.eigenDecompositionHashD.get(num2str([Tf channels]));
%                 else
%                     d_V = eig(self.getNoiseCovarianceMatrix(Tf, channels));
%                     assert(~any(imag(d_V)), 'Covariance matrix has imaginary eigenvalues!');
%                     assert(~any(d_V<0), 'Covariance matrix has negative eigenvalues!');
%                     self.eigenDecompositionHashD.put(num2str([Tf channels]), d_V);
%                 end                
%             else
%                 if self.eigenDecompositionHashV.containsKey(num2str([Tf channels]))
%                     D = diag(self.eigenDecompositionHashD.get(num2str([Tf channels])));
%                     d_V = self.eigenDecompositionHashV.get(num2str([Tf channels]));
%                 else
%                     [d_V, D] = eig(self.getNoiseCovarianceMatrix(Tf, channels));
%                     assert(~any(imag(diag(D))), 'Covariance matrix has imaginary eigenvalues!');                    
%                     assert(~any(diag(D)<0), 'Covariance matrix has negative eigenvalues!');
%                     self.eigenDecompositionHashD.put(num2str([Tf channels]), diag(D));
%                     self.eigenDecompositionHashV.put(num2str([Tf channels]), d_V);
%                 end                
%             end
%         end
        %%% ------------------------------------------------------
        function Cinv = inv(self, varargin)
            C = self.getNoiseCovarianceMatrix(varargin{:});
            Cinv = diag(1./diag(C));
        end
        %%% ------------------------------------------------------
        function ch = sqrt(self, varargin)
            C = self.getNoiseCovarianceMatrix(varargin{:});
            ch = diag(sqrt(diag(C)));
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