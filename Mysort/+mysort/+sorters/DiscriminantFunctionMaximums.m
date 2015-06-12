function [DMax Ids OvpIndex] = DiscriminantFunctionMaximums(D, M, maxTau, resampleP)
    OvpIndex = [];
    ovpRange = -maxTau:maxTau;
    nShifts = length(ovpRange);
    nF = size(D,2);
    L  = size(D,1);
    Tf = floor(size(M,1)/2)+1;   
    ovpIdx = 1;
    DMax = D;
    Ids = zeros(size(D));
    if nargin == 3 || isempty(resampleP) || resampleP == 1
        % build the overlap index
        for f1 = 1:nF
            Ids(:,f1) = f1;
            for f2 = 1:nF
                if f1==f2
                    continue
                end
                for shiftIdx = 1:nShifts        
                    shift = ovpRange(shiftIdx);
                    OvpIndex(ovpIdx,:) = [nF+ovpIdx f1 f2 shift shiftIdx];  
                    shift = ovpRange(shiftIdx);

                    tmp = mysort.util.shiftSubtract(D(:,f1)', -D(:,f2)', shift, 2)' - M(Tf+shift, f1, f2);

                    [DMax(:,f1) maxIds] = max([DMax(:,f1) tmp],[], 2);
                    tmp = Ids(:,f1);
                    tmp(maxIds==2) = OvpIndex(ovpIdx,1);
                    Ids(:,f1) = tmp;
                    ovpIdx = ovpIdx +1;
                end
            end
        end
    else
        ovpRange = resampleP*(-maxTau:maxTau);
        nShifts = length(ovpRange);
        M = mysort.wf.tResample(M, resampleP, 1);
        D = resample(D, resampleP, 1);
        for f1 = 1:nF
            Ids(:,f1) = f1;
            for f2 = 1:nF
                if f1==f2
                    continue
                end
                for shiftIdx = 1:nShifts        
                    shift = ovpRange(shiftIdx);
                    OvpIndex(ovpIdx,:) = [nF+ovpIdx f1 f2 shift shiftIdx];  
                    shift = ovpRange(shiftIdx);

                    tmp = mysort.util.shiftSubtract(D(:,f1)', -D(:,f2)', shift, 2)' - M(Tf*resampleP+shift, f1, f2);
                    tmp = reshape(tmp, resampleP, size(DMax,1))';
                    
                    [DMax(:,f1) maxIds] = max([DMax(:,f1) tmp],[], 2);
                    tmp = Ids(:,f1);
                    tmp(maxIds>1) = OvpIndex(ovpIdx,1);
                    Ids(:,f1) = tmp;
                    ovpIdx = ovpIdx +1;
                end
            end
        end       
        OvpIndex(:,4) = round(OvpIndex(:,4)/resampleP); 
    end