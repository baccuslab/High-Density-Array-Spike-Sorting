function [smad, commonNoiseEpochs, notNoiseEpochsSet] = estimateSigma(X, Tf, thr)
% Compute the median absolute deviation and the corresponding std ersimate
% for every row of data matrix X individually.
% The method calculates an estimate of the noise std by using the MAD. With
% that estimate a thresholding is performed and possible spikes detected.
% On the rest of the signal the MAD computation is repeated.
% Input: 
%    X - data matrix. rows correspond to recording channels
%   Tf - length around a possible spike that will be treated as "not noise"
%  thr - threshold factor for multiplication with the estimated standard
%        deviation of the noise to detect possible spikes.
%        WARNING: Be aware on how to set thr! If you have 100 channels and
%                 you set Tf to 200 and thr to 3.5, the probability for any
%                 given sample to be considered as noise is 
%                 PdetOneSided      = 1 - normcdf(3.5);                
%                 PdetDoubleSided   = 2*PdetOneSided;
%                 PnoDetDoubleSided = 1-PdetDoubleSided;
%                 PnoDetInWindow    = PnoDetDoubleSided^(100*200);
%                 PAtLeastOneDetInWindow = 1-PnoDetInWindow   % = .9999  
% Output:
%              smad - The std as guessed by the MAD per channel
% commonNoiseEpochs - Noise epochs that were common over all channels
% notNoiseEpochsSet - Epochs per channel that were classified as "not noise"

if ~exist('Tf', 'var')
    Tf = 15;
end
maxIter = 5;
if ~exist('thr', 'var')
    thr = 4.5;
end

[L nC] = size(X);

noiseIdx = true(1, L);
smad = zeros(1,nC);

for i=0:maxIter
    if sum(noiseIdx) < .05*L
        warning('Common noise epochs are less than 5% of original data! Taking whole data')
        noiseIdx = true(1, L);
        break
    end
    smad_ = std(X(noiseIdx,:), [], 1);
    smad_(smad_<=10*eps) = 1000;
    if all(smad_==smad | smad_==1000)
        noiseIdx = true(1, L);
        break
    end
    
    smad = smad_;
    noiseIdx = all(abs(X)<repmat(smad,L,1)*thr,2);
    if ~any(noiseIdx)
        warning('No noise epochs left to estimate sigma! Taking whole data');
        noiseIdx = true(1, L);
        break
    end
end
smad = mysort.util.mad(X(noiseIdx,:)')';
smad(smad<=10*eps) = 1000;
if nargout > 1
    commonNoiseEpochs = mysort.epoch.fromBinaryVector(noiseIdx);
    notNoiseEpochsSet = mysort.epoch.removeShort(mysort.epoch.flip(commonNoiseEpochs, L), Tf);
end

