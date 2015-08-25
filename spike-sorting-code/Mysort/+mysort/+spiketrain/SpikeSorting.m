classdef SpikeSorting < handle
    
    properties
        gdf
        unitIDs
        Units
        spikeCounts
        trialTimes
        
        DataSource
    end
    
    methods
        % -----------------------------------------------------------------
        function self = SpikeSorting(gdf, trialTimes, footprints, datasource)
            % Input
            %   gdf   - matrix with 2 columns
            %           first column : ID of neuron
            %           second column: spike time with respect to some global reference time point of 0
            %  trialTimes - (optional) matrix with two columns defining
            %               start- (first column) and end-times of trials
            %  footprints - (optional) foot prints for each unit. 3 matrix
            %               with (time x electrodes x units)
            %  datasource - (optional)
            
            s = size(gdf,2);
            assert(s >= 2 && s <= 3, 'Must have 2 or three columns!');
            if s == 2
                gdf(:,3) = 1;
            end
            if nargin < 2
                trialTimes = [];
            end
            
            self.trialTimes = trialTimes;
            self.gdf = gdf;
            self.unitIDs = unique(gdf(:,1));
            
            if nargin > 2
                assert(size(footprints,3) == length(self.unitIDs), 'There must be one footprint per Unit!');
            else
                footprints = [];
            end
            if nargin > 3
                self.DataSource = datasource;
            else 
                self.DataSource = [];
            end
            
            self.Units = mysort.spiketrain.Unit.empty(0,1);
            st = self.getSpikeTrains();
            for ui=1:length(self.unitIDs)
                U = self.unitIDs(ui);
                if ~isempty(footprints)
                    self.Units(end+1) = mysort.spiketrain.Unit(U, st(:,ui), self, footprints(:,:,ui));
                else
                    self.Units(end+1) = mysort.spiketrain.Unit(U, st(:,ui), self);
                end
            end
            
            self.getSpikeCounts();
            
        end
        
        
        % -----------------------------------------------------------------
        function plot(self)

        end
        
        % -----------------------------------------------------------------
        function gdf = getGdf(self)
            gdf = self.gdf;
        end     
        
        % -----------------------------------------------------------------
        function U = getUnit(self, unitIDs)
            idx = ismember(self.unitIDs, unitIDs);
            U = self.Units(idx);
        end
        
        % -----------------------------------------------------------------
        function st = getSpikeTrains(self, unitIDs)
            if nargin < 2
                unitIDs = self.unitIDs;
            else
                assert(all(ismember(unitIDs, self.unitIDs)), 'Unit not in this sorting!');
            end
            mgdf = self.getMGdf(unitIDs);
            % Remove all spikes that were not in any trial
            mgdf(mgdf(:,3)==0,:) = [];
            
            st = mysort.spiketrain.mgdf2spikeTrains(mgdf);            
        end     
        % -----------------------------------------------------------------
        function spikeCounts = getSpikeCounts(self)
            if isempty(self.spikeCounts)
                [Uq Nq Gq] = unique(self.gdf(:,1)); % Uq: unique rows, Uq = II(Nq,:), II = Uq(Gq,:)
                [self.spikeCounts x] = hist(Gq, length(Uq));
            end
            spikeCounts = self.spikeCounts;
        end
        
        % -----------------------------------------------------------------
        function T = getT_merged(self, units_idx)
            T = [];
            for uidx = units_idx
                T(:,:,uidx) = self.Units(uidx).footprint;
            end
        end
        
        % -----------------------------------------------------------------
        function mgdf = getMGdf(self, unitIDs)
            if isempty(self.trialTimes)
                mgdf = self.getGdf;
                mgdf(:,3) = 1;
                return
            end
            gdf_ = self.getGdf;
            if nargin == 2
                % remove units that were not requested
                gdf_(~ismember(gdf_(:,1), unitIDs), :) = [];
            end
            mgdf = mysort.spiketrain.gdf2multiSessionGdf(gdf_, self.trialTimes);
        end            
    end
end