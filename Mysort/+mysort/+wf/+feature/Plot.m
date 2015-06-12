classdef Plot < guiutil.AmplitudeSelectableAxes
    properties
       % GUI members
      
       
       % internal members
       featureList
    end
    
    methods
        %------------------------------------------------------------------
        function self = Plot(featureList, varargin)
            P.dummy = 1;
            P.callback = []; 
            [P uP] = mysort.util.parseInputs(P, varargin, 'split');
            uP = mysort.util.deflateP(uP);
            self = self@guiutil.AmplitudeSelectableAxes(P.callback, uP{:}, 'clickableTicks', 1:length(featureList));
            
            self.featureList = featureList;
        end
        %------------------------------------------------------------------
        function n = getNFeatures(self)
            n = length(self.featureList);
        end
        %------------------------------------------------------------------
        function update(self)
            nF = self.getNFeatures();
            for i=1:nF
                f = self.featureList(i);
                x = f.get();
%                 amps = G.WF.getFeatures('MinMax');
%                 amps = -amps(:,1);
%                 idx = amps(:,1) > thr;
%                 if ~any(idx)
%                     return
%                 end
%                 chans = G.WF.eventChans(idx);
%                 % plot brushable plots and link data
%                 h = plot(ax, chans, amps(idx), ...
%                     'kx', 'HitTest', 'off', 'markersize',4);
%                 maxi = max(amps(idx));
%                 set(ax, 'xlim', [0 G.WF.nC+1]);
%                 set(ax, 'ylim', [0 maxi*1.01]);
%                 xlabel(ax, 'channel index');
%                 ylabel(ax, 'amplitude [std noise]');
%                 set(ax, 'color', 'none');
            end
        end        
    end
end