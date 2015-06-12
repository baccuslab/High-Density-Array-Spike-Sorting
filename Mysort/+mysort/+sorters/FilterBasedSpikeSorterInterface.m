classdef FilterBasedSpikeSorterInterface < mysort.sorters.OnlineSpikeSorterInterface ...
                                         & mysort.sorters.GaussianNoiseSpikeSorterInterface ...
                                         & mysort.ds.DataSourceInterface
    properties 
        T
        F
        CONF
        Y
        nF        
    end
    
    methods(Abstract)
        fout = getFilterOutputs_(X, filterindex)
    end
    
    methods
        %%% ----------------CONSTRUCTOR---------------------------     
        function self = FilterBasedSpikeSorterInterface(Covest, Tf, T, varargin)
            self = self@mysort.sorters.OnlineSpikeSorterInterface(varargin{:});            
            self = self@mysort.sorters.GaussianNoiseSpikeSorterInterface(Covest, Tf, varargin{:});
            nC = size(Covest.CCol, 2);
            ME = mysort.ds.MultiElectrode((1:nC)', (1:nC)');
            self = self@mysort.ds.DataSourceInterface('FilterBasedSpikeSorterInterface', [], ME);
            mysort.util.checkForForbiddenParameter(varargin, {'Tf'});
            self.P.dummy = [];
            self.P = mysort.util.parseInputs(self.P, varargin);
            
            self.Y = [];
            self.setTemplates(T, nC);
            self.F = [];
            self.nF = [];            
            self.CONF = [];            
        end
      
        %%% ------------------------------------------------------
        function fOut = getData_(self, timeindex, filterindex)
%             X = self.DH(timeindex(1):timeindex(end),:)';
%             fout = self.getFilterOutputs_(X, filterindex)';
            t1 = min(timeindex);
            t2 = max(timeindex);
            [sorting fOut] = self.sortData(t1, t2);
            fOut = fOut(:,timeindex-t1+1)';
        end
        
        %%% ------------------------------------------------------
%         function varargout = size(self, varargin)
%             dims = [size(self.DH,1) self.nF];
%             varargout = matlabfilecentral.parseSize.parseSize(dims, nargout, varargin{:}); 
%         end
        %%% ------------------------------------------------------
        function n = getNSamples_(self)
            n = size(self.DH,1);
        end
      
        %%% ------------------------------------------------------
        function setTemplates(self, T, nC)
            if nargin == 2; nC = 1; end
            assert(ndims(T)==2, ['Templates must be a matrix containing'...
                   ' Templates as rows. If templates have multiple '...
                   'channels, they need to be concatinated']);
            Tf = size(T,2)/nC;
            assert(Tf == round(size(T,2)/nC), ...
                   ['Templates must be a matrix containing'...
                   ' Templates as rows. If templates have multiple '...
                   'channels, they need to be concatinated']);               
%            self.setChannelNumber(nC)
            self.setTemplateLength(Tf)
            self.nF = size(T,1);
            self.T = T;
            self.bReadyToSort = 1;
        end
        
        %%% ------------------------------------------------------
        function setTemplateLength(self, Tf)
            if ~isempty(self.Tf)
                assert(self.Tf==Tf, 'Template dimension is inconsistent');
            else
                self.Tf = Tf;
            end
        end 
        
        
        %%% ------------------------------------------------------
        function plotTemplates(self)
            mysort.plot.spikes(self.T, 'classes', 1:self.nF, 'nC', size(self.DH,2));
            mysort.plot.figureName('Templates');
        end
        
        %%% ------------------------------------------------------
        function RES = plotConfusionMatrix(self)
            mysort.plot.figure('color','w');%'ymax', 1000, 'ymin', -1000,
            RES = mysort.plot.XIvsF(mysort.wf.v2t(self.T,size(self.DH,2)), mysort.wf.v2t(self.F,size(self.DH,2)),...
                             'figure', 0, ...
                            'XIvsF',self.CONF,'title',0,'axistight',0, 'nC', size(self.DH,2));    
            mysort.plot.figureName('ConfusionMatrix');                        
        end        
              
    end
end