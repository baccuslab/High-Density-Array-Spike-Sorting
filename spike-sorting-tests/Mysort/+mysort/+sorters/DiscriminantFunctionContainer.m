classdef DiscriminantFunctionContainer < handle
    properties (SetAccess=private)
    end
    properties
        D
        M
        maxTau
        ovpRange
        DOvp
        DOvpIndex
        DOvpIndexHash
    end
    
    methods
        %%% ----------------CONSTRUCTOR---------------------------
        function self = DiscriminantFunctionContainer(D, M, maxTau)
            self.D = D;
            self.M = M;
            self.maxTau = maxTau;
            self.ovpRange = -maxTau:maxTau;
            nShifts = length(self.ovpRange);
            nF = size(D,2);
            L  = size(D,1);
            nOvpFunctions = nF*(nF-1)/2;
            self.DOvp = zeros(L, nOvpFunctions*nShifts);
            self.DOvpIndex = zeros(nOvpFunctions*nShifts, 5);
            Tf = floor(size(M,1)/2)+1;   
            ovpIdx = 1;
            for f1 = 1:nF
                for f2 = 1:nF
                    if f1==f2
                        continue
                    end
                    for shiftIdx = 1:nShifts
                        shift = self.ovpRange(shiftIdx);
                        self.DOvpIndex(ovpIdx,:) = [nF+ovpIdx f1 f2 shift shiftIdx];
                        tmp = mysort.util.shiftSubtract(D(:,f1)', -D(:,f2)', shift, 2)' - M(Tf+shift, f1, f2);
                        self.DOvp(:,ovpIdx) = tmp; % THIS IS MUCH FASTER USING tmp !! Why? (3min instead of 40min on Harris data d11221.002)
                        % TODO: DEPENDS ON HOW M is defined!!! Might be
                        % M(Tf-shift, f1, f2) or, equivalently, M(Tf+shift,
                        % f2, f1)
                        % Checked, and, no, +shift is correct!
                        ovpIdx = ovpIdx +1;
                    end
                end
            end
            self.DOvpIndexHash = self.TripleHash(self.DOvpIndex(:,2), self.DOvpIndex(:,3), self.DOvpIndex(:,5));
        end
        %--------------------------------------------------------
        function d = get(self, f1, f2, shift)
            [a shiftIdx] = ismember(shift, self.ovpRange);
            hash = self.TripleHash(f1, f2, shiftIdx);
            d = self.DOvp(:,self.DOvpIndexHash == hash);
        end
       
%     end
%     methods(Static)
        %--------------------------------------------------------
        function hash = TripleHash(self, f1, f2, shiftIdx)
%             if f2 < f1
%                 t = f1;
%                 f1 = f2;
%                 f2 = t;
%                 shiftIdx = 
%             end
            hash = f1(:)*10^6 +f2(:)*10^3 + shiftIdx(:);
        end
    end
end
    