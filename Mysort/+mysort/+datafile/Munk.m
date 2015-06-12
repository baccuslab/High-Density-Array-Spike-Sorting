
classdef Munk < mysort.datafile.DataFileInterface
    properties
        channels
        dataFiles
        channelIdx2Id
        channelId2Idx
    end
    
    methods
        %%% ------------------------------------------------------
        function self = Munk(dirName, varargin)
            self = self@mysort.datafile.DataFileInterface(dirName, varargin{:});
        end
        %%% ------------------------------------------------------
        function [bLoaded nC Len] = init(self, dirName)
            % fname is a directory containing the singe channel files
            fnames = dir(dirName);
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
    end
end