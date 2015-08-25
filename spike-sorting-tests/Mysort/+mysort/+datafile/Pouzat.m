
classdef Pouzat < mysort.datafile.DataFileInterface
    properties
        channels
        reference
    end
    
    methods
        %%% ------------------------------------------------------
        function self = Pouzat(dirName, varargin)
            self = self@mysort.datafile.DataFileInterface(dirName, varargin{:});
        end
        %%% ------------------------------------------------------
        function [bLoaded nC Len] = init(self, dirName)
            self.reference = [];
            % Check if data is already converted
            [bLoaded nC Len] = self.init_(dirName);
        end
        %%% ------------------------------------------------------
        function [bLoaded nC Len] = init_(self, fname)
            % fname is a directory containing the singe channel files
            fnames = dir(fname);
            self.channels = [];
            cIdx = 1;
            for i=1:length(fnames)
                if fnames(i).isdir || ...
                   isempty(strfind(fnames(i).name, 'PK')) || ...
                   ~strcmp(fnames(i).name(end-2:end), '.fl') 
                    continue
                end
                
                if ~isempty(strfind(fnames(i).name, 'ref'))
                    self.debugout(sprintf('Found Reference file: %s', fnames(i).name), self.LEVEL_PROCESSING);
                    self.reference = [fname filesep fnames(i).name];
                else
                    self.debugout(sprintf('Found file: %s', fnames(i).name), self.LEVEL_PROCESSING);

                    self.channels(cIdx).filename =  fnames(i).name;
                    tok = regexp(fnames(i).name, 'PK_([0-9]*)\.fl', 'tokens');
                    self.channels(cIdx).Id = str2double(tok{1});
                    self.channels(cIdx).fullfile = [fname filesep fnames(i).name];
                    %channels(i).h5info = hdf5info([fname filesep fnames(i).name]);
                    memmap = memmapfile(self.channels(cIdx).fullfile, 'Format', 'single');
                    self.channels(cIdx).Len = length(memmap.data);
                    cIdx = cIdx +1;
                end
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
            
            Len = self.channels(1).Len;
        end
        
        %%% ------------------------------------------------------
        function X = getData_(self, start, stopp)
            nC = self.nC;
            X = zeros(nC, stopp-start+1);
            cIdx = 1;
            for i = 1:length(self.channels)
                myIdx = i;
                myId  = i;
                memmap = memmapfile(self.channels(cIdx).fullfile, 'Format', 'single');
                X(cIdx,:) = memmap.data(start:stopp);
                cIdx = cIdx +1;
            end
        end 
        
        %%% ------------------------------------------------------
        function ref = getReference(self)
            if ~isempty(self.reference)
                memmap = memmapfile(self.reference, 'Format', 'single');
                ref = memmap.data();    
            else
                error('No reference file was found during initialization!');
            end
        end
        %%% ------------------------------------------------------
        function [locs, pks, x] = getReferenceSpikeTrain(self)
            x = self.getReference();
            [pks, locs] = findpeaks(-x+.001*randn(size(x)),'minpeakheight',.2,...
                                       'minpeakdistance', 45);
            pks = -pks;
        end        
    end
end