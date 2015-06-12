
classdef mteo < mysort.detectors.DetectorInterface
    % Spike detection via Multi resolution Teager Energy Operator (MTEO).
    %
    % out = mcMteoSpikeDetection(X,p)
    % 
    % OUTPUT
    %
    %     out.spikeTimeIndices
    %     out.mteoSignal 
    %     out.smoothedMteoSignal 
    %     out.threshold         The actually used threshold (threshFactor*std(mteo))
    %     out.parameters        The parameters used in the computation. See optional Parameters
    %
    % PARAMETERS
    % 
    %     X       data matrix of size [nrChannels timeSteps]
    %     
    % OPTIONAL PARAMETERS
    %     p.k                     A vector containing the indices of the k-TEOs
    %                             that shall be used in this MTEO, default is k
    %                             = [1 3 5]
    %     p.threshFactor          Determines the threshold of spike detection.
    %                             Treshold will be the std of the nced signal
    %                             times threshFactor.
    methods
        %%% ----------------CONSTRUCTOR---------------------------
        function self = mteo(varargin)
            self = self@mysort.detectors.DetectorInterface(varargin{:});
            self.P.k = [1 3 5];
            self.P.threshFactor = 4;
            self.P = mysort.util.parseInputs(self.P, 'mteo', varargin);    
        end
        %%% ------------------------------------------------------
        function energy = computeEnergy(self, X)
            K = self.P.k;
            verbose = self.P.verbose;
            [nC L] = size(X);

            % compute the k-teo for each k value
            energy = zeros(length(K),L);
            for i=1:length(K)
                k = K(i);

                if verbose
                    disp(['mteo: ' num2str(k) '-teo']);
                end
                teo = mysort.detectors.kteo('k', k);
                sig = teo.detect(X);
                if nC > 1
                    energy(i,:) = max(sig);
                else
                    energy(i,:) = sig;
                end
            end

            % if there were more than one k, take the max over the responses
            if length(K) > 1
                energy = max(energy);
            end            
        end
    end
end




