classdef GuiManager < handle
    properties
        % data
        WF
        T
        
        % gui handles
        normHistPlotThresholdHandle
        
        % iteractives
        normHistThresh
        selTemplateId
        normHistIdx
        
        % vars
        chanAmpThresh
        chanAmpChanIdx
        
        % events
        featuresChanged
        brushingChanged
        selectionChanged
        cutSpikesChanged
        selectedIdxChanged
        templateChanged
        spikesAligned
        templateSelectionChanged 
        chanAmpThreshWasChanged   
        templateWasAccepted
        normHistThreshWasChanged
        normHistSelectionChanged
    end
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = GuiManager(WF, T, varargin)
            self.WF = WF;
            self.T = T;
            ids = T.getTemplateIds;
            self.selTemplateId = ids(1);
            self.normHistPlotThresholdHandle = [];
            self.chanAmpThresh = [];
            self.chanAmpChanIdx = [];            
            self.setUpdated();
        end
        %------------------------------------------------------------------
        function setUpdated(self)
            self.chanAmpThreshWasChanged = false; 
            self.featuresChanged = 0;
            self.brushingChanged = 0;
            self.selectionChanged = 0;
            self.cutSpikesChanged = 0;
            self.selectedIdxChanged = 0;
            self.templateChanged = 0;
            self.spikesAligned = 0;                   
            self.templateSelectionChanged = 0;
            self.templateWasAccepted = 0;
            self.normHistThreshWasChanged = 0;
            self.normHistSelectionChanged = 0;
        end        
        %------------------------------------------------------------------
        function b = spikesChanged(self)
           b =  self.brushingChanged  || self.selectionChanged || ...
                self.cutSpikesChanged || self.selectedIdxChanged || ...
                self.templateChanged  || self.spikesAligned;            
        end
        
        %------------------------------------------------------------------
        function acceptCurrentTemplate(self)
            if self.currentTemplateIsAcceptable()
                st1 = self.WF.eventTimes;
                st2 = self.getCurrentSpikeTrain();
                M = mysort.spiketrain.findNearSpikes(st1, st2, 15);
                if ~isempty(M)
                    f = self.WF.getFeature('MinMax');
                    f.F(M(:,1), 1) = -1;
                end
                self.selTemplateId = self.T.accept(self.selTemplateId);
                self.templateSelectionChanged = 1;
                self.templateWasAccepted = 1;
            end
        end
        %------------------------------------------------------------------
        function evts = getCurrentSpikeTrain(self)
            t = self.getCurrentTemplate();
            idx = t.getIdx();
            evts = self.T.wfs.eventTimes(idx);
        end
        
        %------------------------------------------------------------------
        function b = currentTemplateIsAcceptable(self)
            t = self.getCurrentTemplate();
            b = ~t.isAccepted() && self.currentTemplateHasSourceWfs();
        end
        %------------------------------------------------------------------
        function b = currentTemplateIsAccepted(self)
            t = self.getCurrentTemplate();
            b = t.isAccepted();
        end        
        %------------------------------------------------------------------
        function b = currentTemplateHasSourceWfs(self)
            b = self.T.hasSourceWfs(self.selTemplateId);
        end
        %------------------------------------------------------------------
        function t = getCurrentTemplate(self)
            t = self.T.getTemplateWithId(self.selTemplateId);
        end
        %------------------------------------------------------------------
        function bidx = getCurrentTemplateSourceIdx(self)
            t = self.getCurrentTemplate();
            bidx = t.sourceWfIdx;
        end
        %------------------------------------------------------------------
        function P = getCurrentTemplateProjections(self)
            t = self.getCurrentTemplate();
            f = self.WF.getFeature(['proj_T' num2str(t.id)]);
            P = f.get(t.getIdx());
        end        
        %------------------------------------------------------------------
        function P = getCurrentTemplateSourceProjections(self)
            t = self.getCurrentTemplate();
            f = self.WF.getFeature(['proj_T' num2str(t.id)]);
            P = f.get(t.getSourceIdx);
        end  
        %------------------------------------------------------------------
        function P = getCurrentTemplateExcludedProjections(self)
            t = self.getCurrentTemplate();
            f = self.WF.getFeature(['proj_T' num2str(t.id)]);
            P = f.get(t.getExcludedIdx);
        end          
        %------------------------------------------------------------------
        function setTemplateSelection(self, tidx)
            tid_old = self.selTemplateId;
            tid_new = self.T.getTemplateIds(tidx);
            if length(tid_old) ~= length(tid_new) || any(tid_old ~= tid_new)
                fprintf('Template selection changed.\n');
                self.templateSelectionChanged = 1;
                self.selTemplateId = tid_new;
            end
        end
      
        %------------------------------------------------------------------
        function setChanAmpThr(self, chan, thr)
            if isempty(self.chanAmpThresh) || self.chanAmpThresh~=thr ...
                    || chan ~= self.chanAmpChanIdx
                self.chanAmpThresh = thr;
                self.chanAmpChanIdx = chan;
                self.chanAmpThreshWasChanged = true;
%                 disp('ChannAmpThres was changed');
            end
        end  

        %------------------------------------------------------------------
        function handles = update(self, handles)
            % Retrieve template selection
            if ~self.templateWasAccepted
                tIdx = hdmeagui.gui.getChosenTemplateIdx(handles);
                self.setTemplateSelection(tIdx);
            end

            % Retrieve brushed spikes selection
            if self.templateSelectionChanged || self.currentTemplateIsAccepted()
                self.brushingChanged = 1;
            else
                bIdx = self.getBrushedIdx();
                self.setBrushedSelection(bIdx);
            end            
            if self.featuresChanged || self.templateWasAccepted %|| self.chanAmpThreshWasChanged
                hdmeagui.view.plotChannelAmp(handles);
            end
            if self.brushingChanged
                self.normHistThresh = [];
                self.normHistIdx = [];
            end
            nIdx = self.getNormHistIdx();
            self.setNormHistSelection(nIdx);
            
            hdmeagui.view.updateTemplateList(handles);    
            t = self.selTemplateId;
            if self.templateChanged || self.spikesChanged() && self.T.hasSourceWfs(t) ...
                    || self.normHistSelectionChanged || self.brushingChanged
                C = handles.CONFIG;
                GC = handles.GUI_CONFIG;
                cla(handles.SpikeAxes);
                cla(handles.NormHistAxes);
                cla(handles.ChipAxes);
                cla(handles.DirectionSelectivityAxes);
                cla(handles.ProjIsiAxes);
                cla(handles.IsiAxes);
                
                if self.currentTemplateHasSourceWfs()
                    ax = handles.SpikeAxes;
                    hdmeagui.view.plotCutSpikes(ax, self.T, t, C, GC);
                
                    ax = handles.NormHistAxes;
                    hdmeagui.view.plotNormHist(ax, self, t, C, GC);
                    
                    ax = handles.ProjIsiAxes; 
                    hdmeagui.view.plotAmplVsISI(ax, self, t, C, GC);
%                 hdmeagui.view.plotFischerProj(handles);   
                    
                    ax = handles.IsiAxes;
                    hdmeagui.view.plotISI(ax, self, t, C, GC);
                    
                    ax = handles.ChipAxes; 
                    hdmeagui.view.plotFootprint(ax, self.T, t, C, GC);
                    
                    ax = handles.DirectionSelectivityAxes;
                    hdmeagui.view.plotDirectionSelectivity(ax, self.T, t, C, GC);
                end
            end

            self.setUpdated();          
        end
                
        %------------------------------------------------------------------
        function nIdx = getNormHistIdx(self)
            nIdx = [];
            if isempty(self.normHistThresh)
                return 
            end
            t = self.getCurrentTemplate();
            sIdx = t.getSourceIdx();
            if isempty(sIdx)
                return
            end
            proj = self.WF.getFeatures(['proj_T' num2str(t.id)], sIdx);
            nIdx = proj > self.normHistThresh;
        end
        %------------------------------------------------------------------
        function setNormHistSelection(self, nIdx)
            if length(self.normHistIdx) ~= length(nIdx) || ...
                    any(self.normHistIdx ~= nIdx)
                fprintf('Normhist selection changed.\n');
                self.normHistSelectionChanged = 1;
                self.normHistIdx = nIdx;
                self.T.setTemplateExcludeIdx(self.selTemplateId, ~nIdx);
            end
        end 
        %------------------------------------------------------------------
        function bIdx = getBrushedIdx(self)
            bIdx = [];
            if isempty(self.chanAmpThresh) 
                disp('ChanAmpThres is empty');
                return
            end
            amps = self.WF.getFeatures4Chan('MinMax', self.chanAmpChanIdx);
            amps = -amps(:,1);
            idx = amps(:,1) > self.chanAmpThresh;
            if ~any(idx)
                return
            end            
            bIdx = find(self.WF.eventChans==self.chanAmpChanIdx);
            bIdx = bIdx(idx);
        end
        %------------------------------------------------------------------
        function setBrushedSelection(self, bidx)
            t = self.getCurrentTemplate();
            b = t.getSourceIdx();
            if length(b) ~= length(bidx) || any(bidx ~= b)
                fprintf('Brushing changed.\n');
                self.brushingChanged = 1;
                self.T.setTemplateSourceIdx(self.selTemplateId, bidx);
            end
        end    
        
        %------------------------------------------------------------------
        function handles = alignSpikes(self, handles)
            temp = self.getCurrentTemplate();
            temp.alignSpikes();
            self.spikesAligned = 1;
            handles = self.update(handles);
        end
    end
end