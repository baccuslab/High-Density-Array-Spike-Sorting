classdef ChipElectrodePlot < handle
    properties
        %defaults
        chipSpecification
        chipDimensionsX
        chipDimensionsY
        
        % gui elements
        Parent
        ChipAxes
        
        % plothandles
        m_hBackgroundElectrodes
        m_hAllSessionElectrodes
        m_hSelectedElectrodes
        m_hActiveElectrodes
        
        % callbacks
        
        % variables
        MEA
    end
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = ChipElectrodePlot(parentHandle, varargin)
            self.Parent = parentHandle;
            self.chipSpecification = nr.hidens_get_all_electrodes(2);
            self.chipDimensionsX = [min(self.chipSpecification.x) max(self.chipSpecification.x)];
            self.chipDimensionsY = [min(self.chipSpecification.y) max(self.chipSpecification.y)];           
            self.MEA = [];  
            p = self.Parent;
            self.ChipAxes = axes('Parent', p);
            set(self.ChipAxes, 'Xlim', self.chipDimensionsX,...
                               'Ylim', self.chipDimensionsY,...
                               'XDir','normal',...
                               'YDir','reverse');
        end
        
        %------------------------------------------------------------------
        function setMEA(self, MEA)
            self.MEA = MEA;
            self.update();
        end
        
        %------------------------------------------------------------------
        function update(self)    
            if ~isempty(self.MEA)
                if isempty(self.m_hBackgroundElectrodes) || ~ishandle(self.m_hBackgroundElectrodes)
                    self.m_hBackgroundElectrodes = plot(self.ChipAxes, ...
                            self.chipSpecification.x, ...
                            self.chipSpecification.y, '.', 'markersize', 1, 'color', [.8 .8 .8]);
                end
                activeME = self.MEA.getMultiElectrode();
                active_el = activeME.getElectrodePositions(); 
                
                selectedME = self.MEA.getSelectedSessionsMergedMultiElectrode();
                sel_el = selectedME.getElectrodePositions();
                sel_el_card = selectedME.getElectrodeCardinality();
                
                mergedME = self.MEA.getAllSessionsMergedMultiElectrode();
                allsession_el = mergedME.getElectrodePositions(); 
                allsession_el_card = mergedME.getElectrodeCardinality();
                
                if ishandle(self.m_hActiveElectrodes)
                    delete(self.m_hActiveElectrodes);
                    self.m_hActiveElectrodes = [];
                end
                if ishandle(self.m_hSelectedElectrodes)
                    delete(self.m_hSelectedElectrodes);
                    self.m_hSelectedElectrodes = [];
                end
                if ishandle(self.m_hAllSessionElectrodes)
                    delete(self.m_hAllSessionElectrodes);
                    self.m_hAllSessionElectrodes = [];
                end          
                
                self.m_hAllSessionElectrodes = mysort.mea.plotArray(...
                        allsession_el(:,1) , allsession_el(:,2), 'axesHandle', self.ChipAxes, ...
                        'cla', false, 'color', [.5 .5 .5]);
                
                maxel = self.MEA.getNSessions();
                multiElIdx = allsession_el_card > 1 & allsession_el_card < maxel;
                if any(multiElIdx)
                    mysort.mea.plotArray(...
                        allsession_el(multiElIdx,1) , allsession_el(multiElIdx,2), 'axesHandle', self.ChipAxes, ...
                        'cla', false, 'color', [0 1 0], 'marker', 'x');
                end                
                self.m_hSelectedElectrodes = mysort.mea.plotArray(...
                    sel_el(:,1) , sel_el(:,2), 'axesHandle', self.ChipAxes, ...
                    'cla', false, 'color', [.1 .1 .1]);
                self.m_hActiveElectrodes = mysort.mea.plotArray(...
                    active_el(:,1) , active_el(:,2), 'axesHandle', self.ChipAxes,...
                    'color', 'b', 'cla', false);                
                constantElIdx = allsession_el_card == maxel & maxel > 1;
                if any(constantElIdx)
                    mysort.mea.plotArray(...
                        allsession_el(constantElIdx,1) , allsession_el(constantElIdx,2), 'axesHandle', self.ChipAxes, ...
                        'cla', false, 'color', [1 0 0], 'marker', 'x');
                end

                set(self.ChipAxes, 'Xlim', self.chipDimensionsX,...
                                   'Ylim', self.chipDimensionsY,...
                                   'XDir','normal',...
                                   'YDir','reverse');
            end
        end
    end
end