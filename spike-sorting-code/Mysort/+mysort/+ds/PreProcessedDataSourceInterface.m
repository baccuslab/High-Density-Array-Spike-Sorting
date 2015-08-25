classdef PreProcessedDataSourceInterface < mysort.ds.ExtendedDataSourceInterface
    properties
        bIsPreprocessed
        bDisablePreprocessing
        bCovIsPrecalculated
        precalculatedCovFileName
        preprocessedBufferFilename
        preprocessedBufferH5info
        preprocessedBufferH5path
        spikeSortingFolder
        
        preprocFolder
        preprocSmad
        preprocDet
        preprocCov
    end
    methods(Abstract)
        getData_(self, idx1, idx2)
    end
    
    methods
        %------------------------------------------------------------------
        function self = PreProcessedDataSourceInterface(preprocessedBufferFilename, preprocessedBufferH5path, h5info, spikeSortingFolder, varargin)
            self = self@mysort.ds.ExtendedDataSourceInterface(varargin{:});
%             assert(ischar(prefilterBufferFilename), 'prefilterBufferFilename must be a filename!');
%             assert(~isempty(prefilterBufferFilename), 'prefilterBufferFilename must not be empty!');
            self.preprocessedBufferFilename = preprocessedBufferFilename;
            self.precalculatedCovFileName = [preprocessedBufferFilename(1:end-3) '_cov.h5'];
            self.preprocessedBufferH5path = preprocessedBufferH5path;  
            self.preprocessedBufferH5info = h5info;
            self.spikeSortingFolder = spikeSortingFolder;
            self.bIsPreprocessed = false;
            self.bCovIsPrecalculated = false;
            self.bDisablePreprocessing = false;
            
            p_str = self.preprocessedBufferH5path(2:end);
            idx = strfind(p_str, '/');
            p_str(idx) = '_';
            self.preprocFolder = fullfile(self.preprocessedBufferFilename(1:end-3), p_str);
            self.preprocSmad = fullfile(self.preprocFolder, 'smad.mat');
            self.preprocDet = fullfile(self.preprocFolder, 'det.mat');
            self.preprocCov = fullfile(self.preprocFolder, 'cov.mat');
            if ~exist(self.preprocFolder, 'file')
                mkdir(self.preprocFolder);
            end
            if exist(self.preprocSmad, 'file') && ...
               exist(self.preprocDet, 'file'); 
                self.bIsPreprocessed = 1;
            end
            if exist(self.preprocCov, 'file')
                self.bCovIsPrecalculated = true;
            end
            self.loadSortings();
        end
        %------------------------------------------------------------------            
        function loadSortings(self)
            if isempty(self.spikeSortingFolder)
                self.SpikeSortingContainers = {};
                return
            end
            flist = dir(fullfile(self.spikeSortingFolder, '*.sorting.h5'));
            for i=1:length(flist)
                ffile = fullfile(self.spikeSortingFolder, flist(i).name);
                self.SpikeSortingContainers{end+1} = mysort.mea.PersistentSpikeSortingContainer(ffile, '/sorting');
            end
        end
   
        %------------------------------------------------------------------            
        function preprocess(self)
            if self.bIsPreprocessed
                return
            end            
            if isempty(self.preprocessedBufferFilename)
                error('No preprocess bufferfilename was specified when this object was created. Set a bufferfilename first before trying to preprocess!');
            end
           
            
            if ~exist(self.preprocFolder, 'file')
                mkdir(self.preprocFolder);
            end

            if ~exist(self.preprocSmad, 'file')
                smad = self.noiseStd();
                save(self.preprocSmad, 'smad');
            end            
            
            if ~exist(self.preprocDet, 'file')
                spikeDetection = struct();
                [times pks] = self.detectSpikes();
                save(self.preprocDet, 'times', 'pks');
            end        
            
            if ~exist(self.preprocCov, 'file')
                Cest = self.getCovest();
                CestS = Cest.toStruct();
                save(self.preprocCov, 'CestS');
                self.bCovIsPrecalculated = true;
            end       

        
            self.bIsPreprocessed = true;            
        end
        %------------------------------------------------------------------
        function b = isPreprocessed(self)
            b = self.bIsPreprocessed;
        end
        %------------------------------------------------------------------
        function Cest = getCovest(self, varargin)
            if self.bCovIsPrecalculated && ~self.bDisablePreprocessing
                Cest = self.getCovestFromBufferFile(varargin{:});
            else
                % if not preprocessed, use super method
                Cest = getCovest@mysort.ds.ExtendedDataSourceInterface(self, varargin{:});
                CestS = Cest.toStruct();
                save(self.preprocCov, 'CestS');
                self.bCovIsPrecalculated = true;              
            end
        end             
        %------------------------------------------------------------------
        function Cest = getCovestFromBufferFile(self, varargin)
            if nargin > 1
                warning('Specific parameters for covest are given but data is loaded from buffer file. Parameters are ignored. Delete buffer file if you want other parameters!');
            end
            CestS = load(self.preprocCov, 'CestS');
            Cest = mysort.noise.Covest2(CestS.CestS);
        end               

        %------------------------------------------------------------------
        function smad = noiseStd(self, varargin)
            if exist(self.preprocSmad, 'file') && ~self.bDisablePreprocessing
                smad = self.getNoiseStdFromBufferFile(varargin{:});
            else
                % if not preprocessed, use super method 
                % TODO WHAT HAPPENS IF FIRST CALL DOES NOT QUERY ALL CHANNELS!!
                smad = noiseStd@mysort.ds.ExtendedDataSourceInterface(self, varargin{:});
                if ~self.bDisablePreprocessing
                    save(self.preprocSmad, 'smad');
                end
            end
        end
        %------------------------------------------------------------------
        function smad = getNoiseStdFromBufferFile(self, c, cidx)
            % TODO WHAT HAPPENS IF FIRST CALL DOES NOT QUERY ALL CHANNELS!!
            smad = load(self.preprocSmad);
            smad = smad.smad;
            if nargin > 2
                smad = smad(cidx);
            end
        end
        %------------------------------------------------------------------
        function [times pks] = detectSpikes(self, varargin)
            if exist(self.preprocDet, 'file') && ~self.bDisablePreprocessing
                [times pks] = self.getDetectSpikesFromBufferFile(varargin{:});
            else
                % if not preprocessed, use super method 
                [times pks] = detectSpikes@mysort.ds.ExtendedDataSourceInterface(self, varargin{:});
                if ~self.bDisablePreprocessing
                    save(self.preprocDet, 'times', 'pks');
                end
            end
        end
        %------------------------------------------------------------------
        function [times pks] = getDetectSpikesFromBufferFile(self, varargin)
            d = load(self.preprocDet);
            times = d.times;
            pks = d.pks;
        end        
    end
end