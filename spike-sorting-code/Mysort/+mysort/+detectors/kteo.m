
classdef kteo < mysort.detectors.DetectorInterface

    methods
        %%% ----------------CONSTRUCTOR---------------------------
        function self = kteo(varargin)
            self = self@mysort.detectors.DetectorInterface(varargin{:});
            self.P.k = 5;
            self.P = mysort.util.parseInputs(self.P, 'kteo', varargin);  
        end
        %%% ------------------------------------------------------
        function energy = computeEnergy(self, X)
            [nC L] = size(X);
            k = self.P.k;
            % apply the k-TEO
            % signal = zeros(nC,L);
            % maxLength = L-k;
            % for n=(k+1):maxLength
            %     signal(:,n) = X(:,n).^2 - X(:,n-k).*X(:,n+k);
            % end

            signal = X.^2 - [X(:,1+k:end) zeros(nC,k)] .* [ zeros(nC,k) X(:,1:(end-k))];

            % smooth the result with a modified hamming window 
            h = hamming(4*k + 1);
            normConstant = sqrt( 3*sum(h.^2) + sum(h)^2 );
            h = h./normConstant;

            energy = filter(h,1,signal')';         
        end
    end
end
