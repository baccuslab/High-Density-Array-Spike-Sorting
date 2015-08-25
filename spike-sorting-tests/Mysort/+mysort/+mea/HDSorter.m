classdef HDSorter < handle
    properties
        P
        tempPrefilteredFiles
        tempPrefilteredFullFiles
        M
    end
    
    methods
        % -----------------------------------------------------------------
        function self = HDSorter(P)
            self.init(P);
        end
        
        % -----------------------------------------------------------------
        function init(self, P)
            assert(isfield(P, 'expName'), 'P must have the field: expName');
            assert(isfield(P, 'mapFile'), 'P must have the field: mapFile');
            assert(exist(P.mapFile, 'file')>0, ['Could not find mapFile: "' P.mapFile '"']);
            assert(isfield(P, 'dataFiles'), 'P must have the field: dataFiles');
            assert(iscell(P.dataFiles), 'dataFiles must be a cell array of fully qualified files (path + file name)');
            assert(length(P.dataFiles)>0, 'dataFiles is empty!');
            assert(ischar(P.dataFiles{1}), 'dataFiles is empty!');
            assert(isfield(P, 'sortingOutTempFolder'), 'P must have the field: sortingOutTempFolder');
            assert(isfield(P, 'finalResultOutFolder'), 'P must have the field: finalResultOutFolder');
            assert(isfield(P, 'sortingName'), 'P must have the field: sortingName');
            if ~exist(P.sortingOutTempFolder, 'file')
                mkdir(P.sortingOutTempFolder)
            end
            if ~exist(P.finalResultOutFolder, 'file')
                mkdir(P.finalResultOutFolder)
            end       
            
            bFileNotFound = false;
            for i=1:length(P.dataFiles)
                if ~exist(P.dataFiles{i}, 'file')
                    bFileNotFound = true;
                    fprintf('Could not find data file: %s\n', P.dataFiles{i});
                end
            end
            assert(~bFileNotFound, 'At least one data file could not be found!');
            
            self.M = load(P.mapFile);
            self.P = P;
        end
        
        % -----------------------------------------------------------------
        function prefilterFiles(self)
            self.tempPrefilteredFiles = {};
            flist = self.P.dataFiles;
            outFolder = self.P.sortingOutTempFolder;
            tempPrefilteredFiles = {};
            tempPrefilteredFullFiles = {};
            map = self.M.map;
            for i=1:length(self.P.dataFiles)
                [a,b,c] = fileparts(flist{i});
                fn  = [b c];
                ffn = fullfile(outFolder, fn);
                tempPrefilteredFullFiles{i} = ffn;
                tempPrefilteredFiles{i} = fn;                
                bConvertFile = true;
                if exist(ffn, 'file')
                    try
                        bIsInProcess = hdf5read(ffn, '/bFileIsInProcess');
            %                 bIsInProcess = 1;
                        if bIsInProcess
                            % This means the file is still in process. We assume no
                            % other user/matlab instance is writing so it is
                            % probably a left over from some other run and we
                            % delete it.
                            delete(ffn); 
                        else
                            % The file was already properly processed, we dont need
                            % to do anything.
                            bConvertFile = false;
                        end
                    catch
                        delete(ffn);
                    end
                end  
                if bConvertFile
                    % Only convert file if it was non existent or deleted
%                     self.preprocessFile(, ffn);
                    mysort.mea.preprocessMea1kH5File(flist{i}, map,...
                        'prefilter', 1, 'outFile', ffn, ...
                        'subtractMeanOverAllChannels', 1);
                end                
            end   
            self.tempPrefilteredFiles = tempPrefilteredFiles;
            self.tempPrefilteredFullFiles = tempPrefilteredFullFiles;
        end
        
        % -----------------------------------------------------------------
        function preprocessFile(self, in, out)
%             mysort.mea.preprocessMea1kH5File(in, self.M.map, 'prefilter', 1, 'outFile', out);
        end
        
        % -----------------------------------------------------------------
        function runSorting(self)
            [gdf_merged T_merged localSorting localSortingID] = ana.startHDSorting(...
                self.tempPrefilteredFullFiles, self.P.sortingOutTempFolder);
            save(fullfile(self.P.finalResultOutFolder, [self.P.sortingName '_results.mat']), 'gdf_merged', 'T_merged', 'localSorting', 'localSortingID');    
        end


%         function exportResults(self)
% exportPath = fullfile(outPath, 'Results');
% if ~exist(exportPath, 'file')
%     mkdir(exportPath)
% end
% R = {};
% for i=1:length(h5FileLists)
%     sortOutPath = fullfile(outPath, ['Config' num2str(i)], sortingName);
%     sourceFile = fullfile(sortOutPath, [sortingName '_results.mat'])
%     assert(exist(sourceFile, 'file')>0, 'File not found!');
%     D = load(sourceFile);
%     
%     % get individual files' length in samples
%     compoundMea = mysort.mea.compoundMea(h5FileLists{i}, 'useFilter', 0, 'name', 'PREFILT');
%     L = compoundMea.X.getAllSessionsLength();
%     
%     % break gdf apart
%     start_end_times = [0 cumsum(L)];
%     mgdf = mysort.spiketrain.gdf2multiSessionGdf(D.gdf_merged, start_end_times);
%     for k=1:length(L)
%         gdf = mgdf(mgdf(:,3)==k,[1:2]);
%         R{k, i} = gdf;
%     end
% end
% save(fullfile(depath, [expName '_resultsForMichele.mat']), 'R');
% end       
    end    
end