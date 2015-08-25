
classdef Robin < mysort.datafile.DataFileInterface
    properties
        channels
        channelIdx2Id
        channelId2Idx
    end
    
    methods
        %%% ------------------------------------------------------
        function self = Robin(dirName, varargin)
            self = self@mysort.datafile.DataFileInterface(dirName, varargin{:});
        end
        %%% ------------------------------------------------------
        function [bLoaded nC Len] = init(self, dirName)
            % Check if data is already converted
            convdir = [dirName filesep 'converted'];
            if exist(convdir, 'dir');
                self.debugout('Converted found, loading...', self.LEVEL_FLOW);
            else
                % Not converted, do that
                self.debugout('No converted found, converting...', self.LEVEL_FLOW);
                self.init_(dirName);
                self.convertData(convdir);
            end 
            [bLoaded nC Len] = self.init_(convdir);
        end
        %%% ------------------------------------------------------
        function [bLoaded nC Len] = init_(self, fname)
            % fname is a directory containing the singe channel files
            fnames = dir(fname);
            self.channels = [];
            self.channelIdx2Id = [];
            cIdx = 1;
            for i=1:length(fnames)
                if fnames(i).isdir || ...
                   isempty(strfind(fnames(i).name, 'ch')) || ...
                   ~strcmp(fnames(i).name(end-3:end), '.mat')
                    continue
                end
                self.debugout(sprintf('Found file: %s', fnames(i).name), self.LEVEL_PROCESSING);
                
                self.channels(cIdx).filename =  fnames(i).name;
                tok = regexp(fnames(i).name, 'ch([0-9]*)\.mat', 'tokens');
                self.channels(cIdx).Id = str2double(tok{1});
                self.channels(cIdx).fullfile = [fname filesep fnames(i).name];
                %channels(i).h5info = hdf5info([fname filesep fnames(i).name]);
                memmap = memmapfile(self.channels(cIdx).fullfile, 'Format', 'int16');
                self.channels(cIdx).Len = length(memmap.data);
                self.channelIdx2Id(cIdx) = self.channels(cIdx).Id;
                cIdx = cIdx +1;
            end
            nC = length(self.channels);
            bLoaded = true;
            if ~(nC > 0)
                warning('No Channels loaded!');
                bLoaded = false;
                Len = -1;
                return
            end
            % Check Channel Lens
            for i=2:nC
                bLoaded = bLoaded && (self.channels(i).Len == self.channels(i-1).Len);
            end
            
            % Build Channel ID to Channel Index map
            self.channelId2Idx = nan(max(self.channelIdx2Id),1);
            for idx = 1:nC
                self.channelId2Idx(self.channelIdx2Id(idx)) = idx;
            end
            Len = self.channels(1).Len;
        end
        
        %%% ------------------------------------------------------
        function X = getData_(self, start, stopp)
            nC = self.nC;
            X = zeros(nC, stopp-start+1);
            cIdx = 1;
            for i = 1:length(self.channelId2Idx)
                if isnan(self.channelId2Idx(i))
                    continue
                end
                myIdx = self.channelId2Idx(i);
                myId  = i;
                memmap = memmapfile(self.channels(cIdx).fullfile, 'Format', 'int16');
                X(cIdx,:) = memmap.data(start:stopp);
                cIdx = cIdx +1;
            end
        end
        
        %%% ------------------------------------------------------
        function convertData(self, converted_path)
            path = self.fname;
            if ~exist(converted_path, 'dir')
                mkdir(converted_path);
            end

            for i=1:length(self.channels)
                fpath= self.channels(i).fullfile;
                fprintf('Converting: %s\n', fpath);
                D = load(fpath, 'allwvs');
                fid = fopen([converted_path filesep 'converted_' self.channels(i).filename],'w');
                fwrite(fid, D.allwvs, 'int16');
                fclose(fid);
            end  
        end
    end
end