classdef BOTMInitOptimizer < handle
    properties
        DS
        spikeTrains
        
        cutGdf
%         cutSpikes
        
        templates
        spikeTrainsPerElectrode
        connectionMatrix
        effectiveTemplates
        XCorrContainer
        filterCoefficients
        botmConstants
        
    end
    properties (SetAccess=protected)
        P
    end
    
    methods
        %%% ----------------CONSTRUCTOR---------------------------
        function self = BOTMInitOptimizer(DS, spikeTrains, varargin)
            %    BIO = mysort.mea.BOTMInitOptimizer(DS, spikeTrains, 'foo', bar)
            % Inputs:
            %   DS      - A DataSource object that gives access to the
            %             filtered data. Needed to compute the noise
            %             covariance(sub-) matrices
            %   spikeTrains - spikeTrains for all detected Neurons. Needed
            %                 to exclude spiking periods from noise
            %                 computation
            % (Optional)
            %  'userField'  -  can be freely set by user (must be a native
            %                  type that can be stored in HDF5
            
            P_.userField = [];                   % can be freely specified by user
            P_.filterUsedForPrefiltering = [];   % is needed to implement the prefilter of the data
                                                 % on the hardware exactly the same way as the 
                                                 % filtering was done in software for the initial
                                                 % spike sorting
            P_.cutleft = 10;       % To compute the templates, spikes must be cut from data. This
                                   % parameter sets where the cutting starts with respect to the 
                                   % timepoint in "spikeTrains".
            P_.filterLength = 45;  % Final length of the filters 
            P_.maxNSpikesForTemplateComputation = 1000; % To save time, do not cut all spikes in 
                                                        % spike trains for template calculation
                                                        
            P_.blockSpikesForNoiseComputationThreshold = 10; % in values of DS
            
            P_.useMaxNNoiseSamples = 100000;
                                                        
            P_.Nn_MAX = 500;       % see chMap function for details
            P_.Nc_MAX = 20;        % see chMap function for details
            P_.NCLIMIT = 6;        % see chMap function for details
            P_ = mysort.util.parseInputs(P_, varargin, 'error');
            
            self.DS = DS;
            if size(spikeTrains,2)>1
                spikeTrains = spikeTrains';
            end
            assert(size(spikeTrains,2)==1, 'Invalid SpikeTrain Cell!');
            self.spikeTrains = spikeTrains;
            self.P = P_;
            
            disp('BOTMOptimizer: Computing Templates')
            self.computeTemplates();
            disp('BOTMOptimizer: Computing SpikeTrainsPerElectrode')
            self.computeSpikeTrainsPerElectrode();
            disp('BOTMOptimizer: Initialising NoiseCovarianceMatrix')
            self.XCorrContainer = mysort.noise.XCorrContainer(DS, P_.filterLength-1, ...
                'spikeTrains', self.spikeTrainsPerElectrode, ...
                'useMaxNNoiseSamples', P_.useMaxNNoiseSamples);
            disp('BOTMOptimizer: Computing ConnectionMatrix')
            self.computeConnectionMatrix();
            disp('BOTMOptimizer: Computing Filters and Noise Covariance Matrices')
            self.computeFilters();
        end
        % -----------------------------------------------------------------
        function computeTemplates(self)
            st = self.spikeTrains;
            for i=1:length(st)
                N = length(st{i});
                if N > self.P.maxNSpikesForTemplateComputation
                    idx = randperm(N, self.P.maxNSpikesForTemplateComputation);
                    st{i} = st{i}(idx);
                end
            end
            gdf = mysort.spiketrain.toGdf(st);
            [L, nC] = size(self.DS);
            
            % remove border spikes
            gdf(gdf(:,2)>L-self.P.filterLength,:) = [];
            
            units = unique(gdf(:,1));
            
            
            self.templates = zeros(self.P.filterLength, nC, length(units));
            for uidx = 1:length(units)
                myidx = gdf(:,1) == uidx;
                wfs_ = self.DS.getWaveform(gdf(myidx,2), self.P.cutleft, self.P.filterLength);        
                wfs_(~any(wfs_~=0,2),:) = [];
                self.templates(:,:,uidx) = mysort.wf.v2t(median(wfs_,1), nC);
%                 nSourceSpikesPerTemplateAndChannel(uidx, 1:size(DSFull,2)) = size(wfs_,1); % necessary for saving!
            end
%             self.cutSpikes = wfs_;
            self.cutGdf = gdf;
        end
        % -----------------------------------------------------------------
        function computeSpikeTrainsPerElectrode(self)
            % This function computes for each template on which electrodes
            % it has so much influence as to block its spikes from the
            % noise computation. Then, for each electrode, all spiketrains
            % are merged that would cause visible spikes on it.
            nC = size(self.DS,2);
            self.spikeTrainsPerElectrode = cell(1, nC);
            
            [mi, ma, mi_idx, ma_idx] = mysort.wf.tMinMaxPerTemplate(self.templates);
            M = max(abs(mi), abs(ma));
            B = M > self.P.blockSpikesForNoiseComputationThreshold;
            for i=1:nC
                idx = find(B(:,i));
                self.spikeTrainsPerElectrode{i} = [];
                for k=1:length(idx)
                    self.spikeTrainsPerElectrode{i} = [self.spikeTrainsPerElectrode{i}; self.spikeTrains{idx(k)}];
                end
                self.spikeTrainsPerElectrode{i} = sort(self.spikeTrainsPerElectrode{i});
            end
        end
        
        % -----------------------------------------------------------------
        function computeConnectionMatrix(self)
            T = mysort.wf.t2v(self.templates);
            nT = size(self.templates,3);
            
            % This needs to be done by Jelena
            chm = self.chMap(T, nT, self.P.Nn_MAX, self.P.Nc_MAX, self.P.NCLIMIT);
            self.connectionMatrix = chm;
        end
        % -----------------------------------------------------------------
        function computeFilters(self)
            maxC = max(self.connectionMatrix(:,1));
            nT = size(self.templates,3);
            
            vF = zeros(nT, self.P.filterLength*maxC);
            vFC = cell(1,nT);
            vT = zeros(nT, self.P.filterLength*maxC);
            vTC = cell(1,nT);
            CM = self.connectionMatrix;
            T = self.templates;
            XCCT = self.XCorrContainer;
            BOTMC = zeros(1,nT);
            parfor i=1:nT
                fprintf('%d / %d\n', i, nT);
                myChannels = CM(i, 2:CM(i,1)+1) + 1;
                myT = mysort.wf.t2v(T(:,myChannels,i));
                vTC{i} = myT;
%                 t_chanEmbed = mysort.wf.t2vce(myT);
                t_chanEmbed = mysort.util.embedTime2embedChan(myT, length(myChannels));
                fce = XCCT.invMul(t_chanEmbed, myChannels);
%                 fte = mysort.wf.vce2vte(fce, length(myChannels));
                fte = mysort.util.embedChan2embedTime(fce, length(myChannels));
                vFC{i} = fte;
                BOTMC(i) = fce(:)'*t_chanEmbed(:);
                if 0
                    %%
                    figure
                    plot(mysort.wf.t2v(myT)', 'x-'), title('Tte')
                    figure
                    plot(t_chanEmbed, 'x-'), title('Tce')
                    figure
                    plot(fce, 'x-'), title('Fce')
                    figure
                    plot(fte, 'x-'), title('Fte')
                    
                end
            end
            self.botmConstants = BOTMC;
            for i=1:nT
                myT = vTC{i};
                fte = vFC{i};
                vT(i, 1:length(myT)) = myT;
                vF(i, 1:length(fte)) = fte;
            end
            self.filterCoefficients = vF;
            self.effectiveTemplates = vT;
        end
        
        % -----------------------------------------------------------------
        function chMap = chMap(self, T, Nn, Nn_MAX, Nc_MAX, NCLIMIT)
        % Mail from Jelena 16.1.15:
        % The idea of Nn_MAX and Nc_MAX in the FPGA is to allocate the hardware in the FPGA to process Nn_MAX neurons and Nc_MAX channels per neuron, but to leave the possibility to use a lower number of neurons and/or channels per neuron.
        % So, Nn_MAX and Nc_MAX are the upper limits of the number of neurons and channels per neuron that a particular, programmed, FPGA can handle, and, you don�t need to reprogram the FPGA in order to use it for a different experiment where you need less or equal number of neurons and channels per neuron. In case you use lower number of neurons and channels per neuron, in the channel map variable and the filters coefficients, the unused fields will be set to 0 (so that the length of this variable stays the same in any experiment, defined by Nn_MAX and Nc_MAX).
        %  
        % *NCLIMIT is the current limit on the number of channels per neuron (NCLIMIT<=Nc_MAX).
        % *Nn_MAX = number of templates (but, if lower Nn used in a particular experiment, some templates are simply set to 0)

            % Code from Jelena, chMap.m
            %%%%% determining relevant channels for neurons %%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % T: multi-electrode template
            % F: multi-electrode filter
            % X: input data
            % NCLIMIT: Number of channels per neuron (alternativelly, we could 
            %          use different numbers of channels per neuron, but we can
            %          see if there's any use of that later...)

            L=size(T,2);
            Nc=size(self.DS,2);   
            Lf=L/Nc;         %single-electrode template length


            templates=T;
            neuronMap=zeros(Nn,Nc);
            pwr=zeros(1,Nc);


            for i=1:Nn
                for j=1:Nc
                     pwr(j)=sum(templates(i,(j-1)*Lf+1:j*Lf).^2)/Lf; 
                end

                for k=1:NCLIMIT
                       [tmp_max,maxId]=max(pwr);
                       pwr(maxId)=0;
                       neuronMap(i,maxId)=1;

                end
            end

            %%%%%%%%% from a different file:

            chMap=zeros(Nn_MAX,Nc_MAX+1);

            %(the hardware is configured for some maximum number of neurons and
            %channels, but a different (lower) number can be set during the runtime.

            for i=1:size(neuronMap,1)
                p=2;
                for j=1:size(neuronMap,2)
                    if neuronMap(i,j)==1
                        chMap(i,p)=j-1;
                        p=p+1;
                    end
                end
                chMap(i,1)=p-2;
            end     
        end        
        
        % -----------------------------------------------------------------
        function [F, connM, botmConsts] = getFilterCoefficientsAndConnectionMatrix(self)
            F = self.filterCoefficients;
            connM = self.connectionMatrix;
            botmConsts = self.botmConstants;
        end            
        
        % -----------------------------------------------------------------
        function saveHWSpikeSorterInitFile(self, fileName)
            % Serializes the object to a file
            S = self.toHWSpikeSorterInitStruct();
            save(fileName, 'S', '-v7.3');
        end
 
        % -----------------------------------------------------------------
        function S = toHWSpikeSorterInitStruct(self)
            S.readme = 'This structure contains all necessary information to start a BOTM run. Created by mysort.mea.BOTMInitOptimizer. Can by loaded by calling BIO = mysort.mea.BOTMInitOptimizer(filename) where "filename" is the name of a file that contains this structure.';
            S.createdOn = date();
            S.P = self.P;

            [coefficients, connectionMatrix, botm_constants] = self.getFilterCoefficientsAndConnectionMatrix();
            
            [Tf, nC, nT] = size(self.templates);
            S.filterLength            = Tf;
            S.numberOfNeurons         = nT;
            S.numberOfElectrodes      = nC;
            S.filter.coefficients     = coefficients;
            S.filter.channelMap       = connectionMatrix;
            S.filter.BOTMconstants    = botm_constants;
        end 
    end
    
    methods (Access=protected)
        
    end
end
    


% Graveyard, might be useful again
%             % VERSION 1 ---------------------------------------------------
%             %    BIO = mysort.mea.BOTMInitOptimizer(filename)
%             % Inputs:
%             %   filename - name of a HDF5 file to which a BOTMInitOptimizer
%             %              object was serialized
%             %
%             % VERSION 2 ---------------------------------------------------
% 
%             % VERSION 1 ---------------------------------------------------
%             if nargin == 1
%                 assert(exist(fullT, 'file')==1, 'File Not found!');
%                 self.load(fullT);
%                 return
%             end
%             
%             % VERSION 2 ---------------------------------------------------     
            
            
            
%         % -----------------------------------------------------------------
%         function load(self, fileName)
%             % Loads the object from a file
%             S = load(fileName);
%             self.fromStruct(S);
%         end  




% function S = createBOTMInitStructure(T, F, C, )
%     % This function creates Structure S that contains all necessary
%     % information about the initial spike sorting to start a BOTM run on
%     % new (or the same old) data
%     %
%     % Inputs:
%     %   T    - 3D matrix (time x channels x items) containing the templates
%     %   F    - same for the filters
%     %   C    - time embedded covariance matrix (blocks for channels)
%     %
%     % Outputs:
%     %   S    - structure containing everything
%     
%     
%     [Tf, nC, nT] = size(T);
%     assert(~any(size(T) ~= size(F)), 'T and F must have same size!');
%     assert(size(C,1) == size(C,2), 'Must be square matrix!');
%     assert(size(C,1) == Tf*nC, 'C does not fit T!')
%     