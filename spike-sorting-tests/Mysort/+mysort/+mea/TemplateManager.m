classdef TemplateManager < handle
    properties             
        TemplateList    % stores the Template objects
    end
    
    methods
        %------------------------------------------------------------------
        function self = TemplateManager(TemplateList)
            if isstruct(TemplateList)
                S = TemplateList;
                self.fromStruct(S);
                return
            end
            self.TemplateList = TemplateList;
        end
        %------------------------------------------------------------------
        function fromStruct(self, S)
            self.TemplateList = mysort.mea.Template.empty();
            for i=1:length(S.TemplateList)
                self.TemplateList(i) = mysort.mea.Template(S.TemplateList(i));
            end
        end
        %------------------------------------------------------------------
        function S = toStruct(self)
            for i=1:length(self.TemplateList)
                S.TemplateList(i) = self.TemplateList(i).toStruct();
            end            
            S.date_ = date();
            S.version = 1;
            S.readme = 'This is a struct derived from the class Template. Dont edit if you dont know what you are doing.';
        end        
        %------------------------------------------------------------------
        function I = getIterator(self, idx)
            if nargin == 2
                I = mysort.util.Iterator(self.TemplateList(idx));
            else
                I = mysort.util.Iterator(self.TemplateList);
            end
        end
        
        %------------------------------------------------------------------
        function [wfs cutleft] = getWaveforms4MultiElectrode(self, ME)
            requestedElNumbers = ME.electrodeNumbers;
            nC = length(requestedElNumbers);
            nT = length(self.TemplateList);
            if nT == 0
                wfs = [];
                cutleft = [];
                return
            end
            I = self.getIterator();
            % check that templates all have common cutleft
            t = I.next();
            cutleft = t.cutLeft;
            cutlength = size(t.waveforms,1);
            while I.hasNext()
                t = I.next();
                assert(cutleft == t.cutLeft, 'cutLeft of all templates has to be identical'); 
                assert(cutlength == size(t.waveforms,1), 'cutLength of all templates has to be identical'); 
            end
            
            wfs = zeros(cutlength, nC, nT);
            
            I = self.getIterator();
            while I.hasNext()
                t = I.next();
                wf = t.getWaveform4MultiElectrode(ME);
                wfs(:,:,I.idx) = wf';
            end
        end
        %------------------------------------------------------------------
        function n = getNSourceSpikes4MultiElectrode(self, ME, tidx)
            requestedElNumbers = ME.electrodeNumbers;
            n = [];
            nC = length(requestedElNumbers);
            nT = length(self.TemplateList);
            if nT == 0
                return
            end
            n = zeros(nT, nC);
            if nargin == 3
                I = self.getIterator(tidx);
            else
                I = self.getIterator();
            end
            while I.hasNext()
                t = I.next();
                n(I.idx,:) = t.getNSourceWaveforms4MultiElectrode(ME);
            end
        end        
        %------------------------------------------------------------------
        function deleteTemplateIdx(self, idx)
            self.TemplateList(idx) = [];
        end
    end
end