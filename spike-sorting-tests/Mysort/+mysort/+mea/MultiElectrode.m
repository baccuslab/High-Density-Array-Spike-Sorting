% classdef MultiElectrode < handle
%     properties
%         elNr
%         elPos
%         nC 
%     end
%     
%     methods(Abstract)
%         
%     end
%     
%     methods
%         %------------------------------------------------------------------        
%         function self = MultiElectrode(elPos, elNr)
%             self.elPos = elPos;
%             self.elNr = elNr; 
%             self.nC = length(self.elNr);
%             assert(self.nC == size(self.elPos,1), 'there must be an electrode position for every electrode!');
%             assert(self.nC > 0 , 'there must be an electrode!');
%             assert(size(self.elPos,2) <= 3 , 'electrode positions must have 2 or 3 dimensions!');
%             assert(size(self.elPos,2) >= 2 , 'electrode positions must have 2 or 3 dimensions!');
%         end
%         
%         %------------------------------------------------------------------        
%         function el = getElectrodeList(self)
%             el = [self.elNr self.elPos];
%         end
%         
%         %------------------------------------------------------------------        
%         function idx = elNr2ElIdx(self, elnr)
%             idx = find(self.elNr==elnr);
%         end
%         
%         %------------------------------------------------------------------        
%         function [I R] = getKNearestElectrodes(self, elIdx, k)
%             x = self.elPos(:,1);
%             y = self.elPos(:,2);
%             [I R] = mysort.mea.nearestElectrodes(x, y, elIdx, k);
%         end
%         %------------------------------------------------------------------
%         function [I R] = getElectrodesCloserThan(self, elIdx, dist)
%             x = self.elPos(:,1);
%             y = self.elPos(:,2);
%             [I R] = electrodesCloserThan(x,y,elIdx,dist);
%         end
%         %------------------------------------------------------------------
%         function d = getDistance(elidx1, elidx2)
%             d = norm(self.elPos(elidx1,:) - self.elPos(elidx2,:));
%         end
%     end
% end