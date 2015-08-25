classdef FilterWrapper < handle
    % This class wraps the filter object design process to fullfill
    % two goals: a) be able to build as many copies of a specific filter as
    % you want without the need to redesign the filter coefficients (this
    % is in principle taken care of by the filter design toolbox since the
    % filter design objects seem to be singletons)
    % b) have a way to initialize that object before knowing if it will be
    % used at all. This allows avoiding designing the filter as long as it
    % is not needed
    properties
        % Init Params
        m_strName
        m_nHighPassFilter
        m_nLowPassFilter
        m_nSamplesPerSecond
        m_nFilterOrder
        m_nFilterType
        
        % Filter Object
        m_objHd
        
        % State vars
        m_bIsDesigned       
    end
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = FilterWrapper(hpf, lpf, samplesPerSecond, filterOrder, name, type)
            if nargin == 1
                s = hpf;
                assert(isstruct(s), 'If only one argument is supplied it must be a struct as returned by getStruct()!');
                self.fromStruct(s);
                return
            end
            self.m_nHighPassFilter = hpf;
            self.m_nLowPassFilter = lpf;
            self.m_nSamplesPerSecond = samplesPerSecond;
            self.m_nFilterOrder = filterOrder;
        
            self.m_bIsDesigned = false;
            if nargin > 4
                self.m_strName = name;
            else
                self.m_strName = 'FilterFactory';
            end
            if nargin > 5
                self.m_nFilterType = type;
            else
                self.m_nFilterType = 'butter';
            end
        end
        %------------------------------------------------------------------
        function hd = getFilterCopy(self)
            if ~self.m_bIsDesigned
                fprintf('Designing filter...');
                self.m_objHd = mysort.mea.filter_design(...
                    self.m_nHighPassFilter,...
                    self.m_nLowPassFilter, ...
                    self.m_nSamplesPerSecond, ...
                    self.m_nFilterOrder, ...
                    self.m_nFilterType);
                self.m_bIsDesigned = true;
                fprintf('done.\n'); 
            end
            hd = self.m_objHd.copy();
        end
        %------------------------------------------------------------------
        function str = getName(self)
            str = self.m_strName;
        end
        %------------------------------------------------------------------
        function str = getString(self)
            str = [num2str(self.m_nHighPassFilter) '_'...
                   num2str(self.m_nLowPassFilter)  '_'...
                   num2str(self.m_nFilterOrder) '_'...
                   self.m_nFilterType];
        end
        %------------------------------------------------------------------
        function s = getStruct(self)
            s.hpf = self.m_nHighPassFilter;
            s.lpf = self.m_nLowPassFilter;
            s.filterOrder = self.m_nFilterOrder;
            s.filterType = self.m_nFilterType;
            s.samplesPerSecond = self.m_nSamplesPerSecond;
            s.name = self.m_strName;
        end    
        %------------------------------------------------------------------
        function self = fromStruct(self, s)
            self.m_nHighPassFilter = s.hpf;
            self.m_nLowPassFilter = s.lpf;
            self.m_nFilterOrder = s.filterOrder;
            self.m_nFilterType = s.filterType;
            self.m_nSamplesPerSecond = s.samplesPerSecond;
            self.m_strName = s.name;
            self.m_objHd = [];
            self.m_bIsDesigned = false;
        end         
    end
end