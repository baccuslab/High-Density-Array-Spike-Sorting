function [T tau mMPmask] = tAlignOnMean(T, varargin)
    P.maxIter = 30;
    P.restrictToChannels = [];
    P.restrictToNMaximalValues = [];
    P.restrictMeanToItems = [];
    P.meanMaskMinLen = [];
    P.maxShiftPerIter = 3;
    P.maxShift = [];
    P.useMedian = 0;
    P = mysort.util.parseInputs(P, varargin, 'error');
    [Tf nC N] = size(T);
    mMPmask = [];
    if isempty(P.restrictToChannels)
        V = mysort.wf.t2v(T);
    else
        V = mysort.wf.t2v(T(:,P.restrictToChannels,:));
        nC = length(P.restrictToChannels);
    end
    if isempty(P.restrictMeanToItems)
        itemidx = 1:N;
    else
        itemidx = P.restrictMeanToItems;
    end
    fun = @(x) mean(x,1);
    if P.useMedian
        fun = @(x) median(x,1);
    end
    tau = zeros(size(V,1),1);
    CTRUNC = 1;
%     CDONTTRUCT = 0;
    Vi = V;
    mMPmask = [];
    for iter = 1:P.maxIter
        M = fun(Vi(itemidx,:)); % mean is faster than median!
        if isempty(P.restrictToNMaximalValues) || length(M) <= P.restrictToNMaximalValues
            maxidx = 1:length(M);
        else
            [maxis maxidx] = sort(abs(M));
            maxidx = maxidx(end-P.restrictToNMaximalValues+1:end);
%             idx = mysort.wf.vSubIdx(M, nC, maxIdx);
        end
        shiftRange = -P.maxShiftPerIter:P.maxShiftPerIter;
        MPmask = zeros(size(M));
        MPmask(:,maxidx) = 1;
        if ~isempty(P.meanMaskMinLen)
            mMPmask = mysort.wf.v2m(MPmask, nC);
            % make at least len 3:
            mMPmask2 = [(mMPmask(:,1) | mMPmask(:,2)) mMPmask(:,1:end-2) | mMPmask(:,2:end-1) | mMPmask(:,3:end) (mMPmask(:,end-1) | mMPmask(:,end)) ];
            % make at least len 5:
            mMPmask2 = [(mMPmask2(:,1) | mMPmask2(:,2)) mMPmask2(:,1:end-2) | mMPmask2(:,2:end-1) | mMPmask2(:,3:end) (mMPmask2(:,end-1) | mMPmask2(:,end)) ];
%             figure; plot(mMPmask(1,:)); hold on; plot(mMPmask2(1,:)+.1, 'g');
            for c = 1:size(mMPmask2,1)
                e = mysort.epoch.fromBinaryVector(mMPmask2(c,:));
                eL = mysort.epoch.length(e);
                e(eL < P.meanMaskMinLen+4,:) = [];
%                 assert(~isempty(e), 'Could not find a single mask epoch!');
                be = mysort.epoch.toIdx(e);
                mMPmask(c,:) = 0;
                mMPmask(c,be)= 1;
            end
%             hold on; plot(mMPmask(1,:)+.2, 'r');
            if any(mMPmask)
                MPmask = mysort.wf.m2v(mMPmask);
            end
        end
        MP = M;
        MP(~MPmask) = 0;
        % figure; plot(M); hold on; plot(MP, 'r');
        MP = mysort.wf.vShift(repmat(MP,length(shiftRange),1), nC, shiftRange, CTRUNC);
        
        % project every waveform on all shifted versions of the mean to find
        % best match
        [mx taui] = max(Vi*MP',[],2); % the masking is computanionally bloodily inefficient this way!
        taui = -taui - (P.maxShiftPerIter+1);
        taui = taui - round(median(taui));
        if ~isempty(P.maxShift)
            taui(abs(taui)>P.maxShift) = 0;
        end
        
        % compute next V
        if any(taui~=0)
            tau = tau + taui;  
        else
            break
        end
        
        if 0
            figure
            iidx = find(taui>7 | taui<-7);
            plot(MP', 'k');
            hold on
            plot(Vi(iidx(1),:))
        end
        
        % always use originals for shifting!
        Vi(taui~=0,:) = mysort.wf.vShift(V(taui~=0,:), nC, tau(taui~=0,:), CTRUNC);        
        nShifts = sum(taui~=0);
        fprintf('alignWaveforms: %d shifted\n', nShifts);
    end
    
    T = mysort.wf.tShift(T, tau, CTRUNC);
end