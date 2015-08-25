
classdef threshold < mysort.detectors.DetectorInterface
% Spike detection via simple thresholding of the normalized input signal X.
% The threshold may be given for each channel individually. 
% By default the signal is squared elementwise before applying the
% specified thresholds. 
%
% out = mcThresholdSpikeDetection(X,p)
% 
% OUTPUT
%
%     out.spikeTimeIndices
%     out.thresholdSignal   The max over the thresholded signals of the
%                           individual channels
%     out.smoothedThresholdSignal 
%                           Smoothed version of out.thresholdSignal
%     out.threshold         The actually used threshold (threshFactor*std(signal))
%     out.parameters        The parameters used in the computation. See optional Parameters
%
% PARAMETERS
% 
%     X       data matrix of size [nrChannels timeSteps]
%     
% OPTIONAL PARAMETERS 
%     
%     p.threshFactor          Determines the threshold of spike detection.
%                             Treshold will be the std of the normalized signal
%                             times threshFactor. 
%     p.method                Method must be one of the following
%                                 'square' - Square the input signal
%                                           elementwise before normalization and
%                                           thresholding
%                                 'abs' - Take the absolute value instead of squaring.
%                                 'none' - Do not preprocess the signal,
%                                 simply apply the treshFactor
%     p.refractoryPeriod      The number of time steps that must be
%                             minimally between two potential spikes
%
% Parameters that are not given by the user are assigned default values.
% See out.parameters for the assigned values.
    methods
        %%% ----------------CONSTRUCTOR---------------------------
        function self = threshold(varargin)
            self = self@mysort.detectors.DetectorInterface(varargin{:});
        end
        
        %%% ----------------CONSTRUCTOR INIT---------------------- 
        function init(self, varargin)
            init@mysort.detectors.DetectorInterface(self, varargin{:});
            self.P.method = 'none';
            self.P = mysort.util.parseInputs(self.P, 'threshold', varargin);    
        end
        %%% ------------------------------------------------------
        function energy = computeEnergy(self, X)
            % preprocess the input signal
            % center the input signal around zero mean
            for c=1:size(X,1)
                X(c,:) = X(c,:) - mean(X(c,:));
            end
            switch self.P.method
                case 'none'
                    energy = X;
                case 'minus'
                    energy = -X;                    
                case 'square'
                    energy = X.^2;
                case 'abs'
                    energy = abs(X);
                otherwise
                    error(['mcThresholdSpikeDetection -> unknown "method" parameter: ' self.P.method]);
            end
        end
    end
end