classdef IndividualChannelSpikeDetector < handle
    properties (SetAccess=private)
        
    end
    properties
        P
        bIsTrained
         
        spikesDetectedDown
        pks_down
        spikesDetectedUp
        pks_up
        smad
        dataBuffer
        chunkCounter 
    end
    
    methods
        %%% ----------------CONSTRUCTOR---------------------------
        function self = IndividualChannelSpikeDetector(varargin)
            self.P.minSamplesToTrain = 100000;
            self.P.thresholdFactor = 3.5;
            self.P.minDist = 20;
            self.P = mysort.util.parseInputs(self.P, varargin, 'error');
            self.dataBuffer = {};
            self.chunkCounter = 0;
            self.spikesDetectedDown = {};
            self.pks_down = {};
            self.spikesDetectedUp = {};
            self.pks_up = {};
        end
        
        %%% ------------------------------------------------------
        function processChunk(self, X)
            if self.bIsTrained
                self.processChunk_(X);
            else
                self.addChunkForTraining(X);
            end            
        end
        %%% ------------------------------------------------------
        function addChunkForTraining(self, X)
            self.dataBuffer{end+1,1} = X;
            
            if self.hasEnoughDataForTraining()
                self.train();
            end
        end
        %%% ------------------------------------------------------
        function b = hasEnoughDataForTraining(self)
            L = sum(cellfun(@(x) size(x,1), self.dataBuffer));
            b = L >= self.P.minSamplesToTrain;
        end
        
        %%% ------------------------------------------------------
        function processChunk_(self, X)
            assert(self.bIsTrained, 'Must be trained first!');
            self.chunkCounter = self.chunkCounter+1;
            nC = size(X,2);
            % DETECT SPIKES
            pks_up = cell(1, nC);
            times_up = cell(1, nC);
            pks_down = cell(1, nC);
            times_down = cell(1, nC);            
            smad_ = self.smad;
            thr = self.P.thresholdFactor;
            mDist = self.P.minDist;
            parfor c=1:nC
                [pks_down{1,c}, times_down{1,c}] = findpeaks(-X(:,c), 'MINPEAKHEIGHT', smad_(c)*thr,...
                                                        'MINPEAKDISTANCE', mDist);
%                 [pks_up{1,c}, times_up{1,c}] = findpeaks(X(:,c), 'MINPEAKHEIGHT', smad_(c)*thr,...
%                                                         'MINPEAKDISTANCE', mDist);
            end            
            self.spikesDetectedDown(self.chunkCounter, :) = times_down;
            self.pks_down(self.chunkCounter, :) = pks_down;               
            self.spikesDetectedUp(self.chunkCounter, :) = times_up;
            self.pks_up(self.chunkCounter, :) = pks_up;            
        end
        
        %%% ------------------------------------------------------
        function [stdown pksdown stup pksup] = getSpikeTimesForChunk(self, idx)
            stdown = self.spikesDetectedDown(idx, :);
            pksdown = self.pks_down(idx, :);            
            stup = self.spikesDetectedUp(idx, :);
            pksup = self.pks_up(idx, :);
        end
        
        %%% ------------------------------------------------------
        function train(self)
            % ESTIMATE NOISE STD
            X = cell2mat(self.dataBuffer);
            smad_ = zeros(1, size(X,2));
            
            parfor c=1:size(X,2)
                DS = mysort.ds.Matrix(X(:,c));
                smad_(c) = DS.noiseStd();
            end
            
            self.smad = smad_;
            self.bIsTrained = true;
            for i=1:length(self.dataBuffer)
                self.processChunk_(self.dataBuffer{i});
                self.dataBuffer{i} = [];
            end
            self.dataBuffer = {};
        end 
    end
end